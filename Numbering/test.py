from renumber_pdb_chothia import renumberpdb

# Specify your PDB file path and the desired output path
pdb_file_path = 'trast_H_pos_1_HLC_83922_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb'
output_path = 'testing.pdb'

# Call the function
renumberpdb(pdb_file_path, output_path, chain_a='A', chain_b='B', scheme='chothia')