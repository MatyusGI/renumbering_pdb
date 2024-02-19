from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from Bio import PDB

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

# Replace 'your_pdb_file.pdb' with the path to your PDB file
pdb_file_path = 'trast_H_pos_1_HLC_83922_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb'
chain_a_sequence = extract_sequence('A', pdb_file_path)

# if chain_a_sequence:
#     print(f"Chain A Sequence: {chain_a_sequence}")
# else:
#     print("Chain A not found or no sequence extracted.")

from anarci import anarci

def number_sequence(chain_a_sequence):
    # Make sure chain_a_sequence is not empty
    if chain_a_sequence:
        # Properly handling the output based on the inspected structure
        results = anarci([('A', chain_a_sequence)], scheme="chothia")

        numbered, details = results[0][0], results[1][0]

        return(numbered, details)

        # for i in range( len( numbered ) ): # Iterate over the identified domains (e.g. for an scfv)
        #     numbering = [ (n, a) for n, a in numbered[i][0] if a != '-' ] # Remove gaps if made (imgt scheme)
        #     yield [ (numbering[ri][0], ri+numbered[i][1]) for ri in range( len( numbering ) ) ], details[i]['chain_type'], details[i]


number, detail = number_sequence(chain_a_sequence)

numbers = []
residues = []

for list_item in number:
    for tuple_item in list_item:
        for item in tuple_item:
            # print(item[0])
            residues.append(item[1])
            num = item[0]
            if num[1] == ' ':
                numbers.append(str(num[0]))
            else:
                numbers.append(str(num[0])+num[1])
        break
    
# print(numbers)
# print(residues)


pdb_io = PDB.PDBIO()
pdb_parser = PDB.PDBParser()
structure = pdb_parser.get_structure("structure", pdb_file_path)

# for model in structure:
#     for chain in model:
#         for i, residue in enumerate(chain.get_residues()):
#             print(residue.id)

for model in structure:
    for chain in model:
        if chain.get_id() == 'A':
            for i, residue in enumerate(chain.get_residues()):
                res_id = list(residue.id)
                res_id[1] = numbers[i]
                
                residue.id = tuple(res_id)
                # print(residue.id)

# pdb_io.set_structure(structure)
# pdb_io.save("renumbered_pdb_chain_a.pdb")
                

from biopandas.pdb import PandasPdb
data = PandasPdb().read_pdb(pdb_file_path)
# print(data.df['ATOM']['residue_number'].head())
# print(data.df['ATOM']['chain_id'] == 'A')


for numb in numbers:
    curr_res_num = 1
    for i, residue in enumerate(data.df['ATOM']['residue_number']):
        # print(numb)
        if residue != curr_res_num:
            break
        else:
            data.df['ATOM']['residue_number'][i] = numb

print(data.df['ATOM']['residue_number'].range(20))
