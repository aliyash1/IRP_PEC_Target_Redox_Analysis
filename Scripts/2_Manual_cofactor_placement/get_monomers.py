# get json files to put into alphafold and get monomers
# from a json file with all PECs, get a single PEC

template_file = 'template.csv'
with open('target.csv', 'r') as csv_file:
    template = np.load_csv(template_file, delimiter=',', dtype=str)
    for row  in template:
        name = row[0]
        sequence = row[1]
        json_data = {"{
            "name": f"{name}",
            "modelSeeds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": f"{sequence}",
                        "count": 1
                    }
                },
               
            ]
        },"
   