from vina import Vina
import glob
import os
import sys
import subprocess
from mol2pdb import *

os.mkdir('results_' + sys.argv[1])
fold= 'results_' + sys.argv[1]
def pep_list(n_pep, l_pep, ref_file, path):
    path = path
    lst = ['A','D','E','F','G','H','I','K','L','N','Q','R','S','T','V','W','Y']
    peptides = []
    for i in range(n_pep):
        peptide = []
        for i in range(l_pep):
            peptide.append(random.choice(lst))
        peptides.append(''.join(peptide))

    
    os.mkdir(path + '/minimized')

    for pep in peptides:
        mol=Chem.MolFromSequence(pep)
        AllChem.EmbedMolecule(mol)
        alt_file = path+'/{}.pdb'.format(pep)
        Chem.MolToPDBFile(mol, alt_file)
        filename = point_mutation(ref_file, alt_file, pep, path)
        pdb_files = PDBFile(filename)
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller = Modeller(pdb_files.topology, pdb_files.positions)

        with open((filename), 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)
        file_path = clean_pdb(filename, path + '/{}_Hs.pdb'.format(pep))
        try:
            ffminimization(file_path, path + '/minimized/{}_Hs_output.pdb'.format(pep))
        except ValueError or ZeroDivisionError:
            pass

        


def point_mutation(ref_file, alt_file, seq, path):
    parser = PDBParser(PERMISSIVE = True, QUIET = True)
    ref_structure = parser.get_structure('peptide_structure', ref_file)
    alt_structure = parser.get_structure('mutated_structure', alt_file)
    ref_atoms = []
    alt_atoms = []
    for (ref_model, alt_model) in zip(ref_structure, alt_structure):
        for (ref_chain, alt_chain) in zip(ref_model, alt_model):
            for (ref_res, alt_res) in zip(ref_chain, alt_chain):
                ref_atoms.append(ref_res['CA'])
                alt_atoms.append(alt_res['CA'])
                
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, alt_atoms)
    super_imposer.apply(alt_model.get_atoms())
    print(alt_model.id, super_imposer.rms)
    io=PDBIO()
    io.set_structure(alt_structure)
    filename = path + '/mut_{}.pdb'.format(seq)
    io.save(filename)
    return filename

pep_list(15,6, 'mod_pep.pdb', fold)
   