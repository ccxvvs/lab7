import sys
import pandas as pd
import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, OperateOnResidueSubset, PreventRepackingRLT
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, NeighborhoodResidueSelector, NotResidueSelector
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
import gc
import os

# 1. 
pyrosetta.init("-mute all")

print("Loading structure and calculating baseline score...")
pose = pose_from_pdb("relaxed_0.pdb")
scorefxn = get_fa_scorefxn()
baseline_score = scorefxn(pose)

print("Loading mutation data...")
df = pd.read_csv("mc4r-dms.tsv", sep="\t")
df = df[df['aa'] != '*']

# Translation
aa_1_to_3 = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

# 2. 
task_id = int(sys.argv[1]) if len(sys.argv) > 1 else 1
total_tasks = int(sys.argv[2]) if len(sys.argv) > 2 else 1

chunk_size = len(df) // total_tasks + 1
start_row = (task_id - 1) * chunk_size
end_row = start_row + chunk_size
my_chunk = df.iloc[start_row:end_row]

print(f"Task {task_id}/{total_tasks}: Processing {len(my_chunk)} mutations...")

# 3. 
output_filename = f"predicted_ddg_results_part_{task_id}.csv"
if not os.path.exists(output_filename):
    with open(output_filename, "w") as f:
        f.write("pos,aa,predicted_ddg\n")

# 4. Run mut
count = 0
for index, row in my_chunk.iterrows():
    pos = int(row['pos'])
    mut_aa_1letter = str(row['aa']).upper()
    
    # Skip unrecognized amino acids
    if mut_aa_1letter not in aa_1_to_3:
        continue
        
    mut_aa_3letter = aa_1_to_3[mut_aa_1letter]
    
    pose_index = pose.pdb_info().pdb2pose(chain='A', res=pos)
    if pose_index == 0:
        continue
        
    mut_pose = pose.clone()
    
    # Apply the mutation using the required 3-letter code
    MutateResidue(pose_index, mut_aa_3letter).apply(mut_pose)
    
    focus_selector = ResidueIndexSelector(str(pose_index))
    nbr_selector = NeighborhoodResidueSelector(focus_selector, 6.0, True)
    not_nbr_selector = NotResidueSelector(nbr_selector)
    
    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), not_nbr_selector))
    
    packer = PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    packer.apply(mut_pose)
    
    mut_score = scorefxn(mut_pose)
    ddg = mut_score - baseline_score
    
  
    with open(output_filename, "a") as f:
        f.write(f"{pos},{mut_aa_1letter},{ddg}\n")
        
    
    count += 1
    if count % 100 == 0:
        print(f"Processed {count} / {len(my_chunk)} mutations...")
        gc.collect()

print(f"Finished! Saved to {output_filename}")
