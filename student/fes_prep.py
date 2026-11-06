#!/usr/bin/env python3
"""
fes_prep.py - prepare protein structures with Fe-S clusters for GROMACS.

Two sub-commands:

  prep    PDB files -> one folder per structure with system.top / system.gro (+ tleap and ParmEd logs)
          Needs AmberTools (tleap) and ParmEd, e.g. the AmberTools conda environment.

  report  after energy minimisation (run_em.sh): collect convergence and Fe-S geometry into one CSV.

Examples
  python fes_prep.py prep  inputs/*.pdb  --out systems
  bash   run_em.sh systems                       (GROMACS environment)
  python fes_prep.py report systems              (writes systems/summary.csv)

Settings (force field, water, which parameter set per cluster type, ...) are in config.ini.
Per-cluster exceptions (e.g. a reduced cluster) go in overrides.csv.
"""
import argparse, configparser, csv, glob, itertools, json, math, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARAMS = os.path.join(HERE, 'parameters')

STD_AA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
          'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'CYX', 'CYM', 'ASH', 'GLH', 'LYN', 'ACE', 'NME'}
WATER = {'HOH', 'WAT', 'DOD', 'H2O', 'TIP', 'TIP3', 'SOL'}
CLUSTER_ELEMENTS = {'FE', 'S'}
S_LIGAND_CUTOFF = 2.9     # A, Fe-S(Cys) and Fe-S(cluster)
NO_LIGAND_CUTOFF = 2.6    # A, Fe-N / Fe-O


# ----------------------------------------------------------------------------------------- utilities
def dist(a, b):
    return math.dist(a, b)


def log(msg):
    print(msg, flush=True)


class Atom:
    def __init__(self, name, xyz, element, record):
        self.name, self.xyz, self.element, self.record = name, xyz, element, record


class Residue:
    def __init__(self, chain, seq, icode, resname, record):
        self.chain, self.seq, self.icode, self.resname, self.record = chain, seq, icode, resname, record
        self.atoms = {}       # name -> Atom (insertion ordered)
        self.new_name = None  # residue name written for tleap
        self.ter_after = False

    @property
    def key(self):
        return f'{self.chain}:{self.seq}{self.icode}'.strip()

    def heavy(self):
        return [a for a in self.atoms.values() if a.element != 'H']


def guess_element(name, elem_field):
    e = elem_field.strip().upper()
    if e:
        return e
    n = name.strip().upper().lstrip('0123456789')
    if n.startswith('FE'):
        return 'FE'
    return n[:1]


def read_pdb(path):
    residues, index = [], {}
    for line in open(path, errors='replace'):
        if line.startswith('ENDMDL'):
            break                                   # first model only
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        alt = line[16]
        if alt not in (' ', 'A', '1'):
            continue
        name, resname = line[12:16].strip(), line[17:21].strip()
        chain, seq, icode = line[21], int(line[22:26]), line[26].strip()
        key = (chain, seq, icode, resname)
        if key not in index:
            r = Residue(chain, seq, icode, resname, line[:6].strip())
            index[key] = r
            residues.append(r)
        r = index[key]
        if name in r.atoms:
            continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        r.atoms[name] = Atom(name, xyz, guess_element(name, line[76:78] if len(line) > 77 else ''), line[:6].strip())
    return residues


def read_missing(path, residues):
    """Unmodelled residues (REMARK 465) split into N-terminal / internal / C-terminal per chain, and the number
    of residues with unmodelled atoms (REMARK 470)."""
    missing, natoms = {}, 0
    for line in open(path, errors='replace'):
        if line.startswith('REMARK 465'):
            p = line[10:].split()
            if len(p) >= 3 and len(p[-3]) == 3 and p[-3].isalpha() and p[-3] != 'RES' and re.fullmatch(r'-?\d+[A-Z]?', p[-1]):
                missing.setdefault(p[-2], []).append(int(re.match(r'-?\d+', p[-1]).group()))
        elif line.startswith('REMARK 470') and re.match(r'REMARK 470\s+(\d+\s+)?[A-Z]{3} [A-Za-z0-9]\s*-?\d+', line):
            natoms += 1
    first, last = {}, {}
    for r in residues:
        if r.record == 'ATOM':
            first.setdefault(r.chain, r.seq)
            last[r.chain] = r.seq
    out = {}
    for ch, nums in missing.items():
        n = [x for x in nums if ch in first and x < first[ch]]
        c = [x for x in nums if ch in last and x > last[ch]]
        i = [x for x in nums if x not in n and x not in c]
        out[ch] = dict(n_term=len(n), internal=len(i), c_term=len(c), internal_residues=i)
    return out, natoms


def write_pdb(path, residues):
    out, n = [], 0
    for i, r in enumerate(residues, start=1):
        rec = 'ATOM  ' if r.record == 'ATOM' else 'HETATM'
        for a in r.atoms.values():
            n += 1
            nm = a.name if len(a.name) == 4 else ' ' + a.name
            out.append(f'{rec}{n % 100000:5d} {nm:<4s} {r.new_name:<4s}{r.chain}{i % 10000:4d}    '
                       f'{a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}  1.00  0.00          {a.element.capitalize():>2s}')
        if r.ter_after:
            out.append('TER')
    out.append('END')
    open(path, 'w', newline='\n').write('\n'.join(out) + '\n')


# ----------------------------------------------------------------------------------------- config
def load_config(path):
    cp = configparser.ConfigParser(delimiters=('=',), inline_comment_prefixes=(';', '#'))
    cp.optionxform = str
    cp.read(path)
    return cp


def load_sets():
    sets = {}
    for f in glob.glob(os.path.join(PARAMS, '*', 'set.json')):
        m = json.load(open(f))
        m['dir'] = os.path.dirname(f)
        sets[m['name']] = m
    return sets


def load_overrides(path):
    """(structure, chain, resseq) -> {'set': ..., 'extra_ligands': ...}; chain/resseq may be '*'."""
    out = {}
    if path and os.path.exists(path):
        for row in csv.DictReader(l for l in open(path) if not l.lstrip().startswith('#')):
            if not (row.get('structure') or '').strip():
                continue
            k = (row['structure'].strip(), (row.get('chain') or '*').strip() or '*', (row.get('resseq') or '*').strip() or '*')
            out[k] = {'set': (row.get('set') or '').strip(), 'extra_ligands': (row.get('extra_ligands') or '').strip().lower()}
    return out


def lookup_override(overrides, structure, chain, resseq):
    for k in ((structure, chain, resseq), (structure, chain, '*'), (structure, '*', resseq), (structure, '*', '*')):
        if k in overrides:
            return overrides[k]
    return {}


# ----------------------------------------------------------------------------------------- cleaning
def clean(residues, cfg, report):
    keep_ions = cfg.get('cleanup', 'keep_ions', fallback='').split()
    drop_unknown = cfg.getboolean('cleanup', 'drop_unknown_hetero', fallback=True)
    his_default = cfg.get('cleanup', 'histidine', fallback='HIE')
    kept = []
    for r in residues:
        # drop hydrogens: tleap rebuilds them from the templates
        r.atoms = {n: a for n, a in r.atoms.items() if a.element != 'H' and a.element != 'D'}
        if not r.atoms:
            continue
        if r.resname in WATER:
            report['_removed'].append(r)
            continue
        if r.resname == 'MSE':                                  # selenomethionine -> methionine
            r.resname, r.record = 'MET', 'ATOM'
            if 'SE' in r.atoms:
                a = r.atoms.pop('SE')
                r.atoms['SD'] = Atom('SD', a.xyz, 'S', 'ATOM')
        elems = {a.element for a in r.atoms.values()}
        if 'FE' in elems and elems <= CLUSTER_ELEMENTS:
            r.kind = 'cluster'
        elif r.resname in STD_AA:
            r.kind = 'protein'
            if r.resname == 'HIS':
                r.resname = his_default
        elif r.resname in keep_ions and len(r.atoms) == 1:
            r.kind = 'ion'
        else:
            if not drop_unknown:
                raise RuntimeError(f'unknown residue {r.resname} {r.key}; set drop_unknown_hetero = yes or parametrise it')
            report['warnings'].append(f'removed non-standard residue/ligand {r.resname} {r.key}')
            report['_removed'].append(r)
            continue
        r.new_name = r.resname
        kept.append(r)
    # chain breaks (TER) after gaps and after every non-protein residue
    for a, b in zip(kept, kept[1:]):
        if a.kind != 'protein' or b.kind != 'protein' or a.chain != b.chain:
            a.ter_after = True
            continue
        if 'C' not in a.atoms or 'N' not in b.atoms or dist(a.atoms['C'].xyz, b.atoms['N'].xyz) > 2.0:
            a.ter_after = True
            report['warnings'].append(f'chain break between {a.resname} {a.key} and {b.resname} {b.key} '
                                      f'(both ends become charged termini)')
    if kept:
        kept[-1].ter_after = True
    return kept


# ----------------------------------------------------------------------------------------- clusters
def find_clusters(residues, report):
    protein = [r for r in residues if r.kind == 'protein']
    clusters = []
    for r in residues:
        if r.kind != 'cluster':
            continue
        fe = [a for a in r.atoms.values() if a.element == 'FE']
        s = [a for a in r.atoms.values() if a.element == 'S']
        bonds = {(f.name, x.name) for f in fe for x in s if dist(f.xyz, x.xyz) < S_LIGAND_CUTOFF}
        ligs, extra = {f.name: [] for f in fe}, []
        for p in protein:
            for a in p.heavy():
                for f in fe:
                    d = dist(f.xyz, a.xyz)
                    if p.resname in ('CYS', 'CYX', 'CYM') and a.name == 'SG' and d < S_LIGAND_CUTOFF:
                        ligs[f.name].append(('CYS', p, 'SG'))
                    elif p.resname in ('HIS', 'HID', 'HIE', 'HIP') and a.name in ('ND1', 'NE2') and d < NO_LIGAND_CUTOFF:
                        ligs[f.name].append(('HIS', p, a.name))
                    elif a.element in ('N', 'O', 'S') and d < NO_LIGAND_CUTOFF and a.name not in ('SG',):
                        extra.append(f'{p.resname}{p.seq}:{a.name}->{f.name} ({d:.2f} A)')
        for x in report['_removed']:
            for a in x.heavy():
                for f in fe:
                    if dist(f.xyz, a.xyz) < NO_LIGAND_CUTOFF:
                        report['warnings'].append(f'cluster {r.resname} {r.key}: removed {x.resname} {x.key} atom {a.name} '
                                                  f'was bound to {f.name} ({dist(f.xyz, a.xyz):.2f} A); this Fe is left without that ligand')
        ncys = sum(1 for v in ligs.values() for k, *_ in v if k == 'CYS')
        nhis = sum(1 for v in ligs.values() for k, *_ in v if k == 'HIS')
        core = f'{len(fe)}Fe{len(s)}S' if s else f'{len(fe)}Fe'
        lig = (f'{ncys}Cys' if ncys else '') + (f'{nhis}His' if nhis else '')
        env = f'{core}_{lig}'
        clusters.append(dict(res=r, fe=fe, s=s, bonds=bonds, ligs=ligs, env=env, extra=extra))
    return clusters


def map_cluster(cl, st):
    """Return (atom_name_map crystal->template, [(protein residue, lib resname, donor, template Fe)])."""
    t_fe = [a['name'] for a in st['core_atoms'] if a['element'] == 'Fe']
    t_s = [a['name'] for a in st['core_atoms'] if a['element'] == 'S']
    t_bonds = {tuple(b) if b[0] in t_fe else (b[1], b[0]) for b in st['core_bonds']}
    t_lig = {f: sorted(l['kind'] for l in st['ligands'] if l['fe'] == f) for f in t_fe}
    c_fe, c_s = [a.name for a in cl['fe']], [a.name for a in cl['s']]
    if len(c_fe) != len(t_fe) or len(c_s) != len(t_s):
        raise RuntimeError('cluster composition does not match the parameter set')
    for pf in itertools.permutations(t_fe):
        mfe = dict(zip(c_fe, pf))
        if any(sorted(k for k, *_ in cl['ligs'][c]) != t_lig[mfe[c]] for c in c_fe):
            continue
        for ps in itertools.permutations(t_s):
            ms = dict(zip(c_s, ps))
            if {(mfe[f], ms[x]) for f, x in cl['bonds']} == t_bonds:
                assign = []
                for c in c_fe:
                    for kind in ('CYS', 'HIS'):
                        crys = [(p, donor) for k, p, donor in cl['ligs'][c] if k == kind]
                        temp = [l for l in st['ligands'] if l['fe'] == mfe[c] and l['kind'] == kind]
                        for (p, donor), l in zip(crys, temp):
                            if donor != l['donor']:
                                raise RuntimeError(f'His {p.key} binds through {donor}, the set needs {l["donor"]}')
                            assign.append((p, l['resname'], l['donor'], mfe[c]))
                return {**mfe, **ms}, assign
    raise RuntimeError('no consistent mapping of the cluster onto the parameter set (check Fe-S distances)')


# ----------------------------------------------------------------------------------------- tleap
def stage_set(st, workdir, protein_ff):
    """Copy a parameter set into the work dir; for ff14SB rename the backbone CA type XC -> CX."""
    files = [st['frcmod'], st['core_mol2']] + [l['lib'] for l in st['ligands']]
    for f in dict.fromkeys(files):
        txt = open(os.path.join(st['dir'], f)).read()
        if protein_ff == 'ff14SB':
            txt = txt.replace('"XC"', '"CX"').replace('XC-', 'CX-').replace('-XC', '-CX')
        open(os.path.join(workdir, f), 'w', newline='\n').write(txt)


def leap_input(cfg, used_sets, bonds, report):
    ff = cfg.get('forcefield', 'protein', fallback='ff19SB')
    water = cfg.get('forcefield', 'water', fallback='opc').lower()
    buf = cfg.getfloat('solvation', 'buffer', fallback=10.0)
    box = {'opc': 'OPCBOX', 'tip3p': 'TIP3PBOX'}[water]
    L = [f'source leaprc.protein.{ff}', f'source leaprc.water.{water}']
    if water == 'opc':
        L.append('loadamberparams frcmod.ionslm_126_opc')
    for st in used_sets:
        L.append('addAtomTypes { ' + ' '.join(f'{{ "{t}" "{e}" "sp3" }}' for t, e in st['atom_types'].items()) + ' }')
        L.append(f'loadamberparams {st["frcmod"]}')
        L.append(f'{st["core_resname"]} = loadmol2 {st["core_mol2"]}')
        for lib in dict.fromkeys(l['lib'] for l in st['ligands']):
            L.append(f'loadoff {lib}')
    L.append('mol = loadpdb clean.pdb')
    L += [f'bond mol.{a}.{an} mol.{b}.{bn}' for a, an, b, bn in bonds]
    L += ['saveamberparm mol dry.prmtop dry.inpcrd',
          f'solvatebox mol {box} {buf}',
          'addIonsRand mol Na+ 0', 'addIonsRand mol Cl- 0',
          'saveamberparm mol system.prmtop system.inpcrd', 'savepdb mol system_leap.pdb', 'quit']
    return '\n'.join(L) + '\n'


def run(cmd, cwd, logfile):
    with open(os.path.join(cwd, logfile), 'w') as fh:
        return subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT).returncode


def prep_one(pdb, outroot, cfg, sets, overrides):
    name = os.path.splitext(os.path.basename(pdb))[0]
    work = os.path.join(outroot, name)
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    report = dict(structure=name, input=os.path.abspath(pdb), status='failed', warnings=[], clusters=[], _removed=[])
    try:
        raw = read_pdb(pdb)
        miss, natoms_missing = read_missing(pdb, raw)
        report['missing_residues'] = miss
        report['residues_with_missing_atoms'] = natoms_missing
        for ch, m in miss.items():
            parts = [f'{m["n_term"]} N-terminal'] * bool(m['n_term']) + [f'{m["internal"]} internal'] * bool(m['internal']) \
                + [f'{m["c_term"]} C-terminal'] * bool(m['c_term'])
            report['warnings'].append(f'chain {ch}: unmodelled residues in the PDB ({", ".join(parts)}); '
                                      f'terminal ones are left out, internal gaps become chain breaks')
        if natoms_missing:
            report['warnings'].append(f'{natoms_missing} residues have unmodelled atoms in the PDB; tleap rebuilds them from templates')
        residues = clean(raw, cfg, report)
        allow_all = cfg.getboolean('options', 'treat_extra_ligands_as_unbound', fallback=False)
        clusters = find_clusters(residues, report)
        if not clusters:
            report['warnings'].append('no Fe-S cluster found')
        defaults = dict(cfg.items('parameter_sets')) if cfg.has_section('parameter_sets') else {}
        used, bonds = {}, []
        idx = {id(r): i for i, r in enumerate(residues, start=1)}
        bound_cys = set()
        for cl in clusters:
            r = cl['res']
            ov = lookup_override(overrides, name, r.chain, f'{r.seq}{r.icode}')
            if cl['extra']:
                msg = f'cluster {r.resname} {r.key}: extra protein ligand(s) {", ".join(cl["extra"])}'
                if allow_all or ov.get('extra_ligands') == 'unbound':
                    report['warnings'].append(msg + ' -> kept NON-bonded (only van der Waals/electrostatics with the cluster)')
                else:
                    cl['env'] += '+other'
                    report['warnings'].append(msg + ' -> no parameter set')
            setname = ov.get('set') or defaults.get(cl['env'])
            info = dict(cluster=r.resname, chain=r.chain, resseq=f'{r.seq}{r.icode}', environment=cl['env'], set=setname)
            report['clusters'].append(info)
            if not setname:
                hint = ('an Fe also binds another protein atom; no set exists for that. If that ligand may be left non-bonded, '
                        'add "<structure>,*,*,,unbound" to overrides.csv (or set treat_extra_ligands_as_unbound = yes in config.ini)'
                        if cl['env'].endswith('+other') else
                        'no set exists for this environment yet, or it is missing from [parameter_sets] in config.ini')
                raise RuntimeError(f'no parameter set for {cl["env"]} (cluster {r.resname} {r.key}): {hint}')
            if setname not in sets:
                raise RuntimeError(f'unknown parameter set {setname}')
            st = sets[setname]
            if st['environment'] != cl['env'].split('+')[0]:
                raise RuntimeError(f'set {setname} is for {st["environment"]}, cluster {r.key} is {cl["env"]}')
            amap, assign = map_cluster(cl, st)
            r.atoms = {amap[n]: Atom(amap[n], a.xyz, a.element, a.record) for n, a in r.atoms.items()}
            r.atoms = {a['name']: r.atoms[a['name']] for a in st['core_atoms']}       # template order
            r.new_name = st['core_resname']
            for p, lib, donor, fe in assign:
                if p.new_name not in (p.resname, lib):
                    raise RuntimeError(f'residue {p.key} coordinates two clusters')
                p.new_name = lib
                bound_cys.add(id(p))
                bonds.append((idx[id(p)], donor, idx[id(r)], fe))
            info['mapping'] = {k: v for k, v in amap.items()}
            info['ligands'] = [f'{p.resname}{p.seq}{p.icode}({p.chain})->{lib}:{donor}-{fe}' for p, lib, donor, fe in assign]
            used[setname] = st
        # disulfides among the remaining cysteines
        cys = [r for r in residues if r.kind == 'protein' and r.resname in ('CYS', 'CYX') and id(r) not in bound_cys and 'SG' in r.atoms]
        for a, b in itertools.combinations(cys, 2):
            if dist(a.atoms['SG'].xyz, b.atoms['SG'].xyz) < 2.5:
                a.new_name = b.new_name = 'CYX'
                bonds.append((idx[id(a)], 'SG', idx[id(b)], 'SG'))
                report['warnings'].append(f'disulfide {a.key}-{b.key}')
        for r in cys:
            if r.new_name == 'CYX' and not any(idx[id(r)] in (x[0], x[2]) for x in bonds):
                r.new_name = 'CYS'
        # terminal residues cannot be cluster ligands (no terminal templates)
        for i, r in enumerate(residues):
            if id(r) in bound_cys and ('OXT' in r.atoms or i == 0 or residues[i - 1].ter_after or r.ter_after):
                raise RuntimeError(f'cluster ligand {r.key} is a chain terminus; not supported by the residue libraries')
        write_pdb(os.path.join(work, 'clean.pdb'), residues)
        pff = cfg.get('forcefield', 'protein', fallback='ff19SB')
        for st in used.values():
            stage_set(st, work, pff)
        open(os.path.join(work, 'leap.in'), 'w', newline='\n').write(leap_input(cfg, list(used.values()), bonds, report))
        tleap = cfg.get('programs', 'tleap', fallback='tleap')
        run([tleap, '-f', 'leap.in'], work, 'tleap.out')
        leaplog = open(os.path.join(work, 'leap.log'), errors='replace').read() if os.path.exists(os.path.join(work, 'leap.log')) else ''
        ll = leaplog.splitlines()
        errs = [l.strip() for l in ll if 'Could not find' in l]
        errs += [ll[i + 1].strip() for i, l in enumerate(ll[:-1]) if l.rstrip().endswith('Error!')]
        if not os.path.exists(os.path.join(work, 'system.prmtop')) or errs:
            raise RuntimeError('tleap failed: ' + '; '.join(dict.fromkeys(errs)) if errs else 'tleap failed, see leap.log')
        # convert to GROMACS
        import parmed
        s = parmed.load_file(os.path.join(work, 'system.prmtop'), os.path.join(work, 'system.inpcrd'))
        q = sum(a.charge for a in s.atoms)
        report['natoms'], report['net_charge'] = len(s.atoms), round(q, 4)
        if abs(q) > 0.01:
            report['warnings'].append(f'net charge after neutralisation is {q:.4f}')
        s.save(os.path.join(work, 'system.top'), format='gromacs', overwrite=True)
        s.save(os.path.join(work, 'system.gro'), overwrite=True)
        report['status'] = 'prepared'
    except Exception as e:                                   # noqa: BLE001 - report every failure per structure
        report['error'] = str(e)
    report.pop('_removed', None)
    json.dump(report, open(os.path.join(work, 'prep_report.json'), 'w', newline='\n'), indent=2, default=str)
    msg = f'{name:20s} {report["status"]:9s} ' + ', '.join(f'{c["cluster"]}{c["resseq"]}{c["chain"]}:{c["environment"]}->{c.get("set")}' for c in report['clusters'])
    if 'error' in report:
        msg += f'  ERROR: {report["error"]}'
    log(msg)
    return report


# ----------------------------------------------------------------------------------------- report
def read_gro(path):
    lines = open(path).read().splitlines()
    n = int(lines[1])
    atoms = []
    for l in lines[2:2 + n]:
        atoms.append((int(l[0:5]), l[5:10].strip(), l[10:15].strip(), tuple(10 * float(x) for x in (l[20:28], l[28:36], l[36:44]))))
    return atoms


def em_result(work):
    out = {}
    logf = os.path.join(work, 'em.log')
    if not os.path.exists(logf):
        return out
    txt = open(logf, errors='replace').read()
    m = re.search(r'Potential Energy\s*=\s*([-\d.eE+]+)', txt)
    f = re.search(r'Maximum force\s*=\s*([-\d.eE+]+)\s+on atom\s+(\d+)', txt)
    s = re.search(r'(converged to Fmax < \S+ in \d+ steps|did not converge[^\n]*|reached the maximum number[^\n]*|'
                  r'Energy minimization has stopped[^\n]*)', txt)
    if m:
        out['epot_kJmol'] = float(m.group(1))
    if f:
        out['fmax'], out['fmax_atom'] = float(f.group(1)), int(f.group(2))
    if s:
        out['em_status'] = s.group(1)
    return out


def geometry(gro_atoms, rep, sets):
    """min/max Fe-S(cluster) distance per cluster after EM, using the set templates."""
    residues, prev = [], None
    for resnr, resn, an, xyz in gro_atoms:              # residue numbers wrap at 100000: segment by order
        if (resnr, resn) != prev:
            residues.append((resn, {}))
            prev = (resnr, resn)
        residues[-1][1][an] = xyz
    out, seen = [], {}
    for c in rep['clusters']:
        st = sets.get(c.get('set'))
        if not st:
            continue
        hits = [atoms for resn, atoms in residues if resn == st['core_resname']]
        k = seen.get(st['core_resname'], 0)
        seen[st['core_resname']] = k + 1
        if k >= len(hits) or not st['core_bonds']:
            continue
        res = hits[k]
        d = [dist(res[a], res[b]) for a, b in st['core_bonds']]
        out.append(dict(cluster=f'{c["cluster"]}{c["resseq"]}{c["chain"]}', fe_s_min=min(d), fe_s_max=max(d)))
    return out


def report_all(root):
    sets = load_sets()
    rows = []
    for f in sorted(glob.glob(os.path.join(root, '*', 'prep_report.json'))):
        rep = json.load(open(f))
        work = os.path.dirname(f)
        row = dict(structure=rep['structure'], prep=rep['status'], error=rep.get('error', ''),
                   clusters=';'.join(f'{c["cluster"]}{c["resseq"]}{c["chain"]}:{c["environment"]}:{c.get("set")}' for c in rep['clusters']),
                   natoms=rep.get('natoms', ''), net_charge=rep.get('net_charge', ''), warnings=len(rep['warnings']),
                   missing_internal=sum(m['internal'] for m in rep.get('missing_residues', {}).values()),
                   missing_termini=sum(m['n_term'] + m['c_term'] for m in rep.get('missing_residues', {}).values()),
                   residues_missing_atoms=rep.get('residues_with_missing_atoms', ''))
        row.update(em_result(work))
        gro = os.path.join(work, 'em.gro')
        if os.path.exists(gro):
            geo = geometry(read_gro(gro), rep, sets)
            if geo:
                row['fe_s_max_after_em'] = round(max(g['fe_s_max'] for g in geo), 3)
                row['fe_s_min_after_em'] = round(min(g['fe_s_min'] for g in geo), 3)
        rows.append(row)
    keys = ['structure', 'prep', 'error', 'clusters', 'natoms', 'net_charge', 'warnings', 'missing_internal',
            'missing_termini', 'residues_missing_atoms', 'em_status', 'epot_kJmol',
            'fmax', 'fmax_atom', 'fe_s_min_after_em', 'fe_s_max_after_em']
    out = os.path.join(root, 'summary.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    log(f'wrote {out} ({len(rows)} structures)')


# ----------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prep', help='prepare PDB files for GROMACS')
    p.add_argument('pdb', nargs='+')
    p.add_argument('--out', default='systems')
    p.add_argument('--config', default=os.path.join(HERE, 'config.ini'))
    p.add_argument('--overrides', default=os.path.join(HERE, 'overrides.csv'))
    r = sub.add_parser('report', help='collect EM results into summary.csv')
    r.add_argument('root', nargs='?', default='systems')
    s = sub.add_parser('sets', help='list available parameter sets')
    a = ap.parse_args()
    if a.cmd == 'sets':
        for n, m in sorted(load_sets().items()):
            log(f'{n:42s} {m["environment"]:16s} {m["charge_state"]:3s} total {m["total_charge"]:+.0f}   {m["source"]}')
        return
    if a.cmd == 'report':
        report_all(a.root)
        return
    cfg = load_config(a.config)
    sets = load_sets()
    overrides = load_overrides(a.overrides)
    os.makedirs(a.out, exist_ok=True)
    reps = [prep_one(f, a.out, cfg, sets, overrides) for f in a.pdb]
    ok = sum(r['status'] == 'prepared' for r in reps)
    log(f'\n{ok}/{len(reps)} structures prepared in {a.out}/')


if __name__ == '__main__':
    main()
