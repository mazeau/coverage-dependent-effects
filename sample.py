from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.build import fcc111, add_adsorbate

# building the slab
slab = fcc111('Ni', size=(2,2,3), vacuum=10.0)

dyn = BFGS()
