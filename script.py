from vina import Vina
import glob
import os
import sys
import subprocess

#os.mkdir('results_' + sys.argv[1])
v = Vina(sf_name='vina')

v.set_receptor('mod_prot.pdbqt')

files = glob.glob('ligs/*.pdb')


for i,file in enumerate(files):
    n_file = file[:-4]+'.pdbqt'
    subprocess.run(['obabel -ipdb {} -opdbqt > {}'.format(file, n_file)], shell=True)
    v.set_ligand_from_file(n_file)
    v.compute_vina_maps(center=[-2, 19, 3], box_size=[30, 30, 30])
    print('here', file)

    energy = v.score()
    print('Score before minimization: %.3f (kcal/mol)' % energy[0])
    

    energy_minimized = v.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose('vina_output/prot_ligand_minimized_{}.pdbqt'.format(n_file[-12:-6]), overwrite=True)

    v.dock(exhaustiveness=8, n_poses=1)
    v.write_poses('vina_output/vina_out_{}.pdbqt'.format(n_file[-12:-6]), n_poses=1, overwrite=True)

#def crossover(a, b):
 #  dihedrals =  get_dihedrals(a)
 #  point = random.randint(0, len(dihedrals)-1)
   