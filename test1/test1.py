from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.build import fcc111, add_adsorbate

### building the slab
lc = 3.912 # lattice constant from doi:10.1103/PhysRev.25.753
# for other surfaces, may need to calc lattice constant
slab = fcc111('Pt', size=(2,2,3), a=lc, vacuum=10.0, periodic=True)

### adding the adsorbate
d = 1.10 # found in https://wiki.fysik.dtu.dk/ase/ase/atoms.html
molecule = Atoms('CO', positions=[(0.,0.,0.),(0.,0.,d)])

### setting calculators
calc = Vasp(xc='PBE', setups='recommended')
slab.set_calculator(calc)
molecule.set_calculator(calc)

### calculate total energies for the systems
#e_slab = slab.get_potential_energy()
#e_CO = molecule.get_potential_energy()

### structure relaxation
## relax the slab
dyn = BFGS(slab, trajectory='slab.traj', restart='slab.pckl') # writing
# trajectory to slab.traj and hessian to slab.pckl 
dyn.run(fmax=0.05) # fmax is convergence criteria for force on all atoms
# smaller fmax is more accurate
#dyn.replay_trajectory('slab_history.traj') # adjusts hessian from iterations
# if more restarting w more than 1 prev .traj file, smoosh together
#ase gui [part1.traj] [part2.traj] -o slab_history.traj

## relax the adsorbate
dyn = BFGS(molecule, trajectory='adsorbate.traj', restart='adsorbate.pckl')
dyn.run(fmax=0.05)

### calculate total energies for the systems
e_slab = slab.get_potential_energy()
e_CO = molecule.get_potential_energy()

### add adsorbate
h = 1.85 # random number from N2Cu.py, maybe bond length?
mask = [atom.tag > 1 for atom in slab] # fixing 2nd and 3rd layers
slab.set_constraint(FixAtoms(mask=mask))
add_adsorbate(slab, molecule, h, 'ontop') # might not go ontop?

### optimize all
dyn = BFGS(slab, trajectory='system.traj', restart='system.pckl')
dyn.run(fmax=0.05)

print('Adsorption energy:', e_slab + e_CO - slab.get_potential_energy())
