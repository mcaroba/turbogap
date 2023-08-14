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

db0 = read("../db2.xyz", index=":")

db = []

for k in range(0, len(db0)):
    atoms = db0[k]
    medoid = atoms.get_array("clmds_medoid")
    del atoms.arrays["clmds_medoid"]
    for i in range(0, len(atoms)):
        if medoid[i]:
            this_medoid = np.full(len(atoms), False)
            this_medoid[i] = True
            atoms2 = atoms.copy()
            atoms2.set_array("clmds_medoid", this_medoid)
            db.append(atoms2)

for k in range(0, len(db)):
    atoms = db[k].copy()
    medoid = atoms.get_array("clmds_medoid")
    surface = atoms.get_array("surface")
    slice = False
    transparency = np.zeros(len(atoms))
    cm = atoms.get_center_of_mass()
    atoms.positions -= cm
    for i in range(0, len(atoms)):
        if medoid[i] and not surface[i]:
            v = atoms[i].position
            for j in range(0, len(atoms)):
                u = atoms[j].position
                if np.dot(v,u) > np.dot(v,v) and j != i:
                    transparency[j] = 0.9
            direction = atoms.positions[i]
            dist = 30.
        elif medoid[i]:
            direction = atoms.positions[i]
            dist = 25.

    try:
        del atoms.arrays["fix_atoms"]
    except:
        pass
    try:
        del atoms.arrays["clmds_medoid"]
    except:
        pass
    try:
        del atoms.arrays["medoid"]
    except:
        pass
    try:
        del atoms.arrays["surface"]
    except:
        pass
    try:
        del atoms.arrays["clmds_cluster"]
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
        if medoid[i]:
            k1 = 0.7; k2 = 0.3
        else:
            k1 = 1.; k2 = 0.
        if atoms[i].symbol == "Pt":
            colors[i] = k1*np.array([0.8,0.8,0.8]) + k2*np.array([1,0,0])
        elif atoms[i].symbol == "Au":
            colors[i] = k1*np.array([1,1,0]) + k2*np.array([1,0,0])
        elif atoms[i].symbol == "H":
            colors[i] = k1*np.array([1,1,1]) + k2*np.array([1,0,0])

    pipeline.source.data.particles_.create_property("Color", data=colors)
    pipeline.source.data.particles_.create_property("Transparency", data=transparency)

    tachyon = TachyonRenderer(shadows=False, direct_light_intensity=1.1)

    pipeline.source.data.cell_.vis.render_cell = False

    vp = Viewport()
    vp.type = Viewport.Type.Perspective
    vp.camera_dir = -direction
    vp.camera_pos = atoms.get_center_of_mass() + direction + direction / np.dot(direction,direction)**0.5 * dist
    vp.render_image(filename = "sites_png/%i.png" % (k+1), size=(120,120), alpha=False, renderer=tachyon)
    pipeline.remove_from_scene()
