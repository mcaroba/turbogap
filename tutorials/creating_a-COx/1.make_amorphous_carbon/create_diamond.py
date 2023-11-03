from ase.io import write
from ase import Atoms
import numpy as np

a = 3.5
cell = [a, a, a]
positions = [[0, 0, 0],
             [a/4, a/4, a/4],
             [0, a/2, a/2],
             [a/4, 3*a/4, 3*a/4],
             [a/2, 0, a/2],
             [3*a/4, a/4, 3*a/4],
             [a/2, a/2, 0],
             [3*a/4, 3*a/4, a/4]]

atoms = Atoms("C8", positions = positions, cell = cell, pbc=True)
# If you want a bigger or smaller system, adjust these numbers below
atoms *= (5,5,5)

# Randomize the velocities at the beginning
vel = 0.01 * (np.random.sample([len(atoms),3])-0.5)
atoms.set_array("vel", vel)

# Change the density
target_density = 2.0 # g/cm^3

V    = atoms.cell.volume * (1e-10)**3 * 100**3        # cm^3
mass = len(atoms) * 12.01 * 1.6605402e-27 * 1000  # g
rho = mass / V  # g/cm^3
print(f"> Atoms: V = {V}cm^3, mass = {mass}g,  density = {rho} g/cm^3")

target_volume = rho * V / target_density / 100**3 / (1e-10)**3

L = target_volume**(1.0/3.0) 

scaled_positions = atoms.cell.scaled_positions(atoms.positions)
atoms.set_cell([L, L, L], scale_atoms = True)

V    = atoms.cell.volume * (1e-10)**3 * 100**3        # cm^3
mass = len(atoms) * 12.01 * 1.6605402e-27 * 1000  # g
rho = mass / V  # g/cm^3
print(f"> Atoms: V = {V}cm^3, mass = {mass}g,  density = {rho} g/cm^3")


write("diamond.xyz", atoms)
