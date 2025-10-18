
import glob
import pandas as pd
import gemmi
import MDAnalysis as mda
#from src import mdaCIF
from MDAnalysis.analysis.distances import distance_array


def compute_openness(structure_path, segment_index):
    # Load your structure and trajectory (or just structure)
    universe = mda.Universe(structure_path)
    res1 = universe.select_atoms(f"segindex {segment_index} and resid 112 and name CA")
    if res1.n_atoms != 1:
        print(f"Failed to select residue 1 with 'segindex {segment_index} and resid 112 and name CA'")
    res2 = universe.select_atoms(f"segindex {segment_index} and resid 284 and name CA")
    if res2.n_atoms != 1:
        print(f"Failed to select residue 2 with 'segindex {segment_index} and resid 284 and name CA'")
    dist = distance_array(res1.positions, res2.positions)[0][0]
    return dist

structure_features = []

# ref structures
print("Loading reference structures ...")
# the numbering on these is not straight-forward, need to do sequence alignment
OPRM1_pdb = pd.read_csv("data/pdb_structures/OPRM1/OPRM1_pdb_summary_20250531.csv")
ref_structures_path = "data/pdb_structures/OPRM1"
for index, row in OPRM1_pdb.iterrows():
    structure_id = row['Entry ID']
    segment_index = row["Entity ID"] - 1
    structure_index = 0
    print(f"Computing features for {structure_id} segment {segment_index} index {structure_index} ...")
    structure_path = f"data/pdb_structures/OPRM1/{structure_id}.cif"
    structure = gemmi.read_structure(structure_path)
    structure.write_pdb('tmp.pdb')
    structure_features.append({
        "source" : "PDB",
        "structure_id" : structure_id,
        "structure_index" : structure_index,
        "structure_path" : structure_path,
        "openness" : compute_openness("tmp.pdb", segment_index = segment_index)})


source = "boltz2_prediction_20250606/OPRM_HUMAN"
print(f"Loading {source} structures ...")
structure_paths = glob.glob(f"intermediate/{source}/*/boltz_results_*/predictions/*/*_model_*.cif")
for index, structure_path in enumerate(structure_paths):
    structure_id = structure_path.split("/")[-1].split("_")[:-2]
    structure_id = "_".join(structure_id)
    structure_index = structure_path.split("/")[-1].split("_")[-1]
    segment_index = 0
    print(f"Computing features for {structure_id} index {structure_index} segment {segment_index}...")
    structure = gemmi.read_structure(structure_path)
    structure.write_pdb('tmp.pdb')
    structure_features.append({
        "source" : source,
        "structure_id" : structure_id,
        "structure_index" : structure_index,
        "structure_path" : structure_path,
        "openness" : compute_openness("tmp.pdb", segment_index = segment_index)})


source = "boltz2_prediction_extra1_20250606/OPRM_HUMAN"
print(f"Loading {source} structures ...")
structure_paths = glob.glob(f"intermediate/{source}/boltz_results_*/predictions/*/*_model_*.cif")
for index, structure_path in enumerate(structure_paths):
    structure_id = structure_path.split("/")[-1].split("_")[:-2]
    structure_id = "_".join(structure_id)
    structure_index = structure_path.split("/")[-1].split("_")[-1]
    segment_index = 0
    print(f"Computing features for {structure_id} index {structure_index} segment {segment_index} ...")
    structure = gemmi.read_structure(structure_path)
    structure.write_pdb('tmp.pdb')
    structure_features.append({
        "source" : source,
        "structure_id" : structure_id,
        "structure_index" : structure_index,        
        "structure_path" : structure_path,
        "openness" : compute_openness("tmp.pdb", segment_index = segment_index)})
    
# remove temporary pdb file, if it exists    
if os.path.exists("tmp.pdb"): os.remove(file_path)
    
structure_features_df = pd.DataFrame(structure_features)

structure_features_df.to_csv(
    "product/structure_features_20250608.tsv",
    sep = "\t",
    index = False)
