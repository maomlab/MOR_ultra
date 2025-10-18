import os
import re
import glob
import pandas as pd
import tqdm
from Bio.PDB import MMCIFParser, Superimposer, Selection
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.cealign import CEAligner

def calculate_rmsd_cealigner(ref, target):
    # Initialize the CEAligner
    ce_aligner = CEAligner()
    # Perform the alignment
    ce_aligner.set_reference(ref)
    ce_aligner.align(target)
    # Get the RMSD from the alignment
    rmsd = ce_aligner.rms
    return rmsd


def calculate_rmsd(ref_atoms, target_atoms):
    """Calculate RMSD after structural alignment."""
    if len(ref_atoms) != len(target_atoms) or len(ref_atoms) == 0:
        print("Number of atoms in the reference and target structure don't match")
        return None
    
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, target_atoms)
    super_imposer.apply(target_atoms)
    return super_imposer.rms



print("Loading reference structures ...")
OPRM1_pdb = pd.read_csv("data/pdb_structures/OPRM1/OPRM1_pdb_summary_20250531.csv")
ref_structures_path = "data/pdb_structures/OPRM1"

parser = MMCIFParser(QUIET=True)
ref_structures = []
for index, row in OPRM1_pdb.iterrows():
    structure_id = row['Entry ID']
    file_path = f"data/pdb_structures/OPRM1/{structure_id}.cif"
    structure = parser.get_structure(structure_id, file_path)
    chains = [chain for chain in structure[0]]
    for chain_index, chain in enumerate(chains):
        if chain_index + 1 != row["Entity ID"]:
            structure[0].detach_child(chain.get_id())
    ref_structures.append({
        'structure_id' : structure_id,
        'file_path' : file_path,
        'structure' : structure})

ref_structures = pd.DataFrame(ref_structures)
print(f"Found {len(ref_structures)} reference structures")

print("Loading target structures ...")
file_path_pattern = re.compile(
        r'data/'
        r'(?P<source>[^_]+)_'                           # source: everything before first underscore
        r'(?P<type>[^_]+)_'                             # type: predictions
        r'(?P<params>[^_]+)_'                           # params: JackHMMER
        r'(?P<date_code>\d{8})/'                        # date_code: 8-digit date
        r'(?P<ligand>.+?)/'                             # ligand: everything until next slash (non-greedy)
        r'pred\.rank_(?P<index>\d+)\.cif$'              # index: number after 'pred.rank_'
    ) 
target_structures = []
for file_path in glob.glob("data/chai1_*/*/*.cif"):
    match = file_path_pattern.match(file_path)
    if match:
        file_path_parsed = match.groupdict()
    else:
        import pdb
        pdb.set_trace()
        raise ValueError("Path {file_path} does not match expected format.")
    structure = parser.get_structure(file_path, file_path)
    structure[0].detach_child('B')
    values = {
        'file_path' : file_path,
        'structure' : structure}
    values.update(file_path_parsed)
    target_structures.append(values)


target_structures = pd.DataFrame(target_structures)    
print(f"Found {len(target_structures)} target structures")

print("Computing RMSDs ...")
rmsds = []
for ref_index, ref_structure in ref_structures.iterrows():
    print(f"Computing RMSDs to {ref_structure['file_path']}")
    for target_index, target_structure in tqdm.tqdm(target_structures.iterrows()):
        rmsd = calculate_rmsd_cealigner(ref_structure['structure'], target_structure['structure'])
        rmsds.append({
            'ref_path': ref_structure['file_path'],
            'target_path': target_structure['file_path'],
            'rmsd': rmsd})
rmsds = pd.DataFrame(rmsds)

rmsds.to_csv("intermediate/rmsd_2.tsv", sep = "\t", index = False)        
