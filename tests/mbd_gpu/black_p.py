"""Black phosphorus (Cmca, 8 atoms per conventional cell) as extended xyz: a
dense system, with hundreds of atoms in every MBD cluster.
Usage: black_p.py <nx> <ny> <nz> <out.xyz>
"""
import sys

A = 3.3136
B = 10.478
C = 4.3763
Y = 0.10168
Z = 0.08056


def cell_fractions():
    base = [(0.0, Y, Z), (0.0, -Y, -Z), (0.0, 0.5 + Y, 0.5 - Z), (0.0, 0.5 - Y, 0.5 + Z)]
    shifted = [(x + 0.5, y + 0.5, z) for (x, y, z) in base]
    return [(x % 1.0, y % 1.0, z % 1.0) for (x, y, z) in base + shifted]


def supercell(nx, ny, nz):
    atoms = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for (x, y, z) in cell_fractions():
                    atoms.append(((x + i) * A, (y + j) * B, (z + k) * C))
    return atoms


def write_xyz(path, atoms, nx, ny, nz):
    lattice = f"{nx * A} 0.0 0.0 0.0 {ny * B} 0.0 0.0 0.0 {nz * C}"
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f'Lattice="{lattice}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for (x, y, z) in atoms:
            f.write(f"P {x:.8f} {y:.8f} {z:.8f}\n")


def main():
    nx, ny, nz = (int(v) for v in sys.argv[1:4])
    write_xyz(sys.argv[4], supercell(nx, ny, nz), nx, ny, nz)


main()
