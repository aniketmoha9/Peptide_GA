from vina import Vina
import glob
import os
import sys
import subprocess
import pandas as pd
from Bio.PDB import *
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem

def point_mutation(ref_file, alt_file, seq):
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
    print(seq)
    filename = base_path+'/mut_{}.pdb'.format(seq)
    io.save(filename)
    return filename

def convert2pdb(seq,alt_file):
    print(seq)
    mol = Chem.MolFromSequence(seq)
    AllChem.EmbedMolecule(mol)
    Chem.MolToPDBFile(mol, alt_file)

def crossover(parent_a,parent_b):
    offspring = parent_a[0:3] + parent_b[-3:]
    alt_file=base_path+'/{}.pdb'.format(offspring)
    convert2pdb(offspring, alt_file)
    ref_file = 'mod_pep.pdb'
    mut_file = point_mutation(ref_file, alt_file, offspring)
    n_file = base_path+'/{}.pdbqt'.format(offspring)
    subprocess.run(['obabel -ipdb {} -opdbqt > {}'.format(mut_file, n_file)], shell=True)
    docking(n_file,i-1)

def save_energies(path):
    files = glob.glob(path+'/gen{}/vina_output/vina_out*.pdbqt'.format(i-1))
    energies = []
    for file in files:
        f = open(file,"r")
        lines = f.readlines()
        for line in lines:

            if line[0:18]=='REMARK VINA RESULT':
                energies.append(line.split()[3])
    d = {'files': files, 'energies': energies}
    df = pd.DataFrame(data=d)
    df = df.sort_values(by=['energies'], ascending=False)
    df.to_csv(path+'/energies{}.csv'.format(i-1))
    return df.files[0], df.files[1]


def docking(n_file,p):
    v.set_ligand_from_file(n_file)
    v.compute_vina_maps(center=[-2, 19, 3], box_size=[30, 30, 30])
    energy = v.score()
    print('Score before minimization: %.3f (kcal/mol)' % energy[0])
    energy_minimized = v.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose(base_path+'/gen{}/vina_output/prot_ligand_minimized_{}.pdbqt'.format(p, n_file[-12:-6]), overwrite=True)
    v.dock(exhaustiveness=8, n_poses=1)
    v.write_poses(base_path + '/gen{}/vina_output/vina_out_{}.pdbqt'.format(p, n_file[-12:-6]), n_poses=1, overwrite=True)



v = Vina(sf_name='vina')

v.set_receptor('mod_prot.pdbqt') #path of the protein pdbqt which needs to be docked

base_path = 'ga_run{}'.format(sys.argv[1])
while os.path.isdir(base_path):
    base_path = 'ga_run{}'.format(int(sys.argv[1])+1)

os.mkdir(base_path)

i = 1
generation=10

files = glob.glob('ligs/prot_ligand_minimized*.pdbqt')
#os.mkdir(base_path + '/gen{}'.format(i))
#os.mkdir(base_path + '/gen{}/vina_output'.format(i))
#files = glob.glob('vina_output/vina_out*.pdbqt')
#import shutil
#for f in fils:
#   shutil.copy(f, base_path + '/gen{}/vina_output'.format(i))
#a,b = save_energies(base_path)
#print(a,b)
while True:
    print(f'Generation {i}')
    os.mkdir(base_path + '/gen{}'.format(i))
    os.mkdir(base_path + '/gen{}/vina_output'.format(i))

    for file in files:
        docking(file,i)
        #shutil.copy(file, base_path + '/gen{}/vina_output'.format(i))
    if i==generation:
        break
    i += 1
    #save_energies(base_path)
#   #df = pd.read_csv(base_path+'/energies.csv')
#   #print(df.files[0], df.files[1])
    parent_a, parent_b = save_energies(base_path)
    crossover(parent_a[-12:-6], parent_b[-12:-6])
    files = glob.glob('{}/gen{}/vina_output/prot_ligand_minimized*.pdbqt'.format(base_path,i-1))

