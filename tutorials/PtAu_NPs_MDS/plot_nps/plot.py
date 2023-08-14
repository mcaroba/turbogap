import sys
from ase.io import read,write
from ase.visualize import view
import numpy as np
from ovito.io.ase import ase_to_ovito
from ovito.pipeline import StaticSource, Pipeline
from ovito.vis import Viewport, BondsVis, ParticlesVis, TachyonRenderer
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier, ColorCodingModifier
from scipy.special import erf

db0 = read("../db1.xyz", index=":")
cl = []
db = []

for atoms in db0:
    if atoms.info["clmds_medoid"]:
        db.append(atoms)
        cl.append(atoms.info["clmds_cluster"])

dbt = np.array(db, dtype=object)
db = dbt[cl]

for k in range(0, len(db)):
    atoms = db[k].copy()
    try:
        del atoms.arrays["fix_atoms"]
    except:
        pass
    atoms_data = ase_to_ovito(atoms)
    pipeline = Pipeline(source = StaticSource(data = atoms_data))
    pipeline.add_to_scene()

    n_Pt = atoms.symbols.count("Pt")
    n_Au = atoms.symbols.count("Au")
    n_H = atoms.symbols.count("H")

    types = pipeline.source.data.particles.particle_types_

    radius = {"Au": 1.8, "Pt": 1.8, "H": 0.5}
    for n in range(1,4):
        try:
            el = types.type_by_id_(n).name
            types.type_by_id_(n).radius = radius[el]
        except:
            pass

    colors = np.zeros([len(atoms),3])
    for i in range(0, len(atoms)):
        if atoms[i].symbol == "Pt":
            colors[i] = np.array([0.8,0.8,0.8])
        elif atoms[i].symbol == "Au":
            colors[i] = np.array([1,1,0])
        elif atoms[i].symbol == "H":
            colors[i] = np.array([1,1,1])

    pipeline.source.data.particles_.create_property("Color", data=colors)

    tachyon = TachyonRenderer(shadows=False, direct_light_intensity=1.1)

    pipeline.source.data.cell_.vis.render_cell = False

    vp = Viewport()
    vp.type = Viewport.Type.Perspective
    direction = (2, 3, -1)
    vp.camera_dir = direction
    vp.camera_pos = atoms.get_center_of_mass() - direction / np.dot(direction,direction)**0.5 * 30.
    vp.render_image(filename = "np_png/%i.png" % (k+1), size=(120,120), alpha=False, renderer=tachyon)
    pipeline.remove_from_scene()

#   Print progress
    sys.stdout.write('\rProgress:%6.1f%%' % ((k + 1)*100./len(db)) )
    sys.stdout.flush()
