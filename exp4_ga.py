from vina import Vina
import glob
import os
import sys
import subprocess
import shutil
import pandas as pd
from Bio.PDB import *
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from mol2pdb import *

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



    pdb_files = PDBFile(mut_file)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = Modeller(pdb_files.topology, pdb_files.positions)

    with open((mut_file), 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    file_path = clean_pdb(mut_file, base_path + '/{}_Hs.pdb'.format(offspring))
    try:
        ffminimization(file_path, base_path + '/minimized/{}.pdb'.format(offspring))
    except ValueError or ZeroDivisionError:
        pass

    n_file = base_path+'/{}.pdbqt'.format(offspring)
    subprocess.run(['obabel -ipdb {} -opdbqt > {}'.format(base_path + '/minimized/{}.pdb'.format(offspring), n_file)], shell=True)
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
    df = df[:population_size]
    df = df.reset_index()
    df = df.drop(['index'], axis=1)
    df.to_csv(path+'/energies{}.csv'.format(i))
    return df.files


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
os.mkdir(base_path+'/ligs')
subprocess.run(['python3 script2.py {}'.format(base_path)], shell=True)
i = 1
generation=10

files = glob.glob(base_path + '/ligs/prot_ligand_minimized*.pdbqt')
files1 = glob.glob(base_path + '/ligs/vina_out*.pdbqt')

os.mkdir(base_path + '/gen{}'.format(i))
os.mkdir(base_path + '/gen{}/vina_output'.format(i))

for file in files:
    shutil.copy(file, base_path+'/gen{}/vina_output'.format(i))

for file1 in files1:
    shutil.copy(file1, base_path+'/gen{}/vina_output'.format(i))



population_size = len(files)





while True:
    print(f'Generation {i}')
    i+=1
    os.mkdir(base_path + '/gen{}'.format(i))
    os.mkdir(base_path + '/gen{}/vina_output'.format(i))
    pop = save_energies(base_path)
    crossover(pop[0][-12:-6], pop[1][-12:-6])
    pop = save_energies(base_path)
    files = pop

    for file in files:
        file = base_path + '/gen{}/vina_output/prot_ligand_minimized_{}.pdbqt'.format(i-1, file[-12:-6])
        docking(file,i)
        #shutil.copy(file, base_path + '/gen{}/vina_output'.format(i))
    if i==generation:
        break

