from nequip.ase.nequip_calculator import NequIPCalculator
from ase.io import *
from ase.calculators.singlepoint import SinglePointCalculator as SPC
from ase.io.trajectory import TrajectoryWriter
from tqdm import tqdm
chemical_symbol_to_type = {
  "Al":"Al",
  'O': "O",
  "Pt": "Pt",
  "C": "C",
  "H": "H",
}

calc = NequIPCalculator.from_deployed_model('./models/md_aug_final_1.pth', device='cuda', species_to_type_name=chemical_symbol_to_type)     
#files_test = ["al2o3_pt8_hx_c3h6_test_perturb.traj", "al2o3_pt8_hx_test.traj", "al2o3_py8_hx_c3h6_perturb.traj", "new_al2o3_pt8_hx_c3h6_test.traj", "perturb_test.traj"]
files_test = ["target_md_aug_full_test.traj"]
# find all files in the current directory that end with train in the name that also end in .xyz
#"al2o3_pt8_hx_c3h6_test_perturb.traj", "al2o3_pt8_hx_test.traj", "al2o3_py8_hx_c3h6_perturb.traj", "new_al2o3_pt8_hx_c3h6_test.traj", "perturb_test.traj"


for file in files_test:
    nn_pred_file_name = file.split('.')[0] + '_nn.traj'
    traj = read(file,':')
    trajw = TrajectoryWriter(nn_pred_file_name,'a')


    for atoms in tqdm(traj):
        atoms.set_calculator(calc)
        en = atoms.get_potential_energy()
        fn = atoms.get_forces()
        atoms.set_calculator(SPC(atoms=atoms, energy = en, forces = fn))
        trajw.write(atoms)
