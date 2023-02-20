from ase.io import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from ase.geometry import wrap_positions
import argparse

def prep_traj(traj):
    data=[]
    for atoms in traj:
        del atoms.constraints
        atoms.set_pbc([1,1,0])
        atoms.set_calculator()
        cell = atoms.get_cell()
        pos = atoms.get_positions()
        pos2 = wrap_positions(pos, cell)
        atoms.set_positions(pos2)
        data.append(atoms)

    return data

def check_in_train(traj, train_traj):
    data_t = traj.copy()
    train_data_t = train_traj.copy()

    print('Prepping Data')
    train_data = prep_traj(train_data_t)
    data = prep_traj(data_t)
    train, valid=[],[]

    print('Creating train/valid list')
    for i in tqdm(range(len(data))):
        atoms = data[i]
        if atoms in train_data:
            train.append(i)
        else:
            valid.append(i)
   
    return train, valid

def plot_parity(traj, traj_nn, train_perturb, train_perturb_nn,  title=None, save=None):
    
    e_bulk_slab = -.13642542E+04
    e_dft = np.array([(atoms.get_potential_energy())/len(atoms) for atoms in traj])
    e_nn =  np.array([(atoms.get_potential_energy())/len(atoms) for atoms in traj_nn])
    e_perturb = np.array([(atoms.get_potential_energy())/len(atoms) for atoms in train_perturb])
    e_perturb_nn = np.array([(atoms.get_potential_energy())/len(atoms) for atoms in train_perturb_nn])

    f_dft = [atoms.get_forces() for atoms in traj]
    f_nn =  [atoms.get_forces() for atoms in traj_nn]
    f_perturb = [atoms.get_forces() for atoms in train_perturb]
    f_perturb_nn = [atoms.get_forces() for atoms in train_perturb_nn]


    ee = abs(e_dft-e_nn)*1000
    e_rmse = np.sqrt(np.mean(ee**2))
    e_all = np.concatenate((e_dft,e_nn, e_perturb, e_perturb_nn))
    e_min = np.min(np.min(e_all))
    e_max = np.max(np.max(e_all))
    xx = [e_min,e_max]
    x_text = np.mean(xx)-0.05*(e_max-e_min)
    y_text = np.mean(xx)-0.3*(e_max-e_min)

    plt.figure(figsize=(10,7.5))

    plt.subplot(1, 2, 1)
    plt.plot(e_perturb,e_perturb_nn,'o', color='#ff7f0e')
    plt.plot(e_dft,e_nn,'o', color='#1f77b4')

    plt.plot(xx,xx,'--')
    plt.text(x_text,y_text,'RMSE(T) = %.02f meV/atom'% (e_rmse))
    plt.xlabel('$E_{DFT}$ eV/atom')
    plt.ylabel('$E_{NN}$ eV/atom')

    f_dft = np.concatenate(f_dft).ravel()
    f_nn = np.concatenate(f_nn).ravel()
    f_perturb = np.concatenate(f_perturb).ravel()
    f_perturb_nn = np.concatenate(f_perturb_nn).ravel()
    fe = abs(f_dft-f_nn)
    f_rmse = np.sqrt(np.mean(fe**2))

    f_all = np.concatenate((f_dft,f_nn))
    f_min = np.min(np.min(f_all))
    f_max = np.max(np.max(f_all))
    xx = [f_min,f_max]
    x_text = np.mean(xx)-0.05*(f_max-f_min)
    y_text = np.mean(xx)-0.3*(f_max-f_min)
    
    plt.subplot(1, 2, 2)
    plt.plot(f_perturb,f_perturb_nn,'o', color='#ff7f0e')
    plt.plot(f_dft,f_nn,'o', color='#1f77b4')
    plt.plot(xx,xx,'--')
    
    
    
    plt.text(x_text,y_text,'RMSE(T) = %.02f eV/$\AA$'% (f_rmse))
    plt.xlabel('$F_{DFT}$ eV/$\AA$')
    plt.ylabel('$F_{NN}$ eV/$\AA$')
    # add title to overall figure
    if title:
        plt.suptitle(title, fontsize=16)
        
    plt.savefig('parity_perturb.png')
    plt.show()
        
if __name__=="__main__":
    #parser = argparse.ArgumentParser()
    #parser.add_argument('--d', help="Ref DFT trajectory", default=None)
    #parser.add_argument('--n', help="NN trajectory", default=None)
    #parser.add_argument('--t', help="Train trajectory", default=None)

    #args = parser.parse_args()

    #traj = read(args.d,':')
    #traj_nn = read(args.n,':')
    traj = "geng_only.traj"
    traj_nn = "geng_only_nn.traj"
    train_perturb = "perturb_test.traj"
    train_perturb_nn = "perturb_test_nn.traj"
    traj = read(traj,':')
    traj_nn = read(traj_nn,':')
    train_perturb = read(train_perturb,':')
    train_perturb_nn = read(train_perturb_nn,':')

    plot_parity(traj, traj_nn, train_perturb, train_perturb_nn, title = "$Al_{2}O_{3}$ + Hx + $Pt_{8}$ + $C_{3}H_{6}$ Perturb Test")
    #plot_parity(traj, traj_nn, title = "$Al_{2}O_{3}$ + Hx + $Pt_{8}$ + $C_{3}H_{6}$ Perturb Test")
