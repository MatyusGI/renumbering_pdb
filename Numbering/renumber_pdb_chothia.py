from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio import PDB
from abnumber import Chain
from biopandas.pdb import PandasPdb

def extract_sequence(chain_id, pdb_file_path):
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file_path)
    sequence = ''

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[0] == ' ':  # Only standard amino acids
                        sequence += seq1(residue.resname, undef_code='X')
                return sequence
    return None

def update_residue_numbers(pdb_data, chain_id, new_residue_numbers):
    # Filter the DataFrame for chain_id
    chain = pdb_data.df['ATOM'][pdb_data.df['ATOM']['chain_id'] == chain_id]
    
    # Update the residue_number for chain_id
    residue_seen = set()
    current_residue_number = 0
    for index, row in chain.iterrows():
        # Check if the residue number has already been updated
        if row['residue_number'] not in residue_seen:
            residue_seen.add(row['residue_number'])
            if current_residue_number < len(new_residue_numbers):
                # Update all atoms of the current residue
                pdb_data.df['ATOM'].loc[(pdb_data.df['ATOM']['chain_id'] == chain_id) & (pdb_data.df['ATOM']['residue_number'] == row['residue_number']), 'residue_number'] = new_residue_numbers[current_residue_number]
                current_residue_number += 1

def renumberpdb(pdb_file_path, output_path='updated_structure_abnumber.pdb', chain_a='A', chain_b='B', scheme='chothia'):
    chain_a_sequence = extract_sequence(chain_a, pdb_file_path)
    chain_b_sequence = extract_sequence(chain_b, pdb_file_path)

    chain_a_number_abnum = Chain(chain_a_sequence, scheme=scheme)
    chain_b_number_abnum = Chain(chain_b_sequence, scheme=scheme)

    numbered_residues_b = [str(pos)[1:] for pos, aa in chain_b_number_abnum]
    numbered_residues_a = [str(pos)[1:] for pos, aa in chain_a_number_abnum]

    # Load your PDB data
    pdb_data = PandasPdb().read_pdb(pdb_file_path)

    # The new residue numbers provided
    new_residue_numbers_a = numbered_residues_a
    new_residue_numbers_b = numbered_residues_b

    update_residue_numbers(pdb_data, chain_a, new_residue_numbers_a)
    update_residue_numbers(pdb_data, chain_b, new_residue_numbers_b)

    # Save the updated DataFrame back to a new PDB file
    pdb_data.to_pdb(path=output_path, records=None, gz=False, append_newline=True)

if __name__ == "__main__":
    # Replace 'your_pdb_file.pdb' with the path to your PDB file
    pdb_file_path = 'trast_H_pos_1_HLC_83922_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb'
    renumberpdb(pdb_file_path)
