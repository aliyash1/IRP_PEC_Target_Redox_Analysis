import pandas as pd
import glob
import pickle

# Load features and add pH
path = glob.glob("extracted_features/features_r5_final_*.csv")[0]
df = pd.read_csv(path, index_col="cofactor_id")
df["pH"] = 7.0

# Load models