from ase.io import read,write

# Read in the database
db = []
for i in range(1,660+1):
    atoms = read("trajs/trajs_1.1/%i.xyz" % i, index=-1)
    db.append(atoms)

write("db0.xyz", db)
