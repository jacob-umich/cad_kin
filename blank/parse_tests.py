import rigid_kin
import parse_cad
import numpy as np
import matplotlib.pyplot as plt

bar = rigid_kin.Structure('program/rigid_kin/5_bar.json')
bar.gen_constraints()
bar.add_parametric_constraint([
    {"c3":[7],"c2":[6]},
    {"c1":[2],"c0":[2]},
])
geometry = bar.cad('program/rigid_kin/5_bar.wls')
labels = bar.get_all_labels(4)
tree = parse_cad.CadTree(geometry,10,4,labels)
branches = tree.branches
mats = out[2].get_constraint_matrix()
labels = bar.get_dof_labels()
fig,ax = plt.subplots(dpi=720)
im = ax.imshow(mat)
cbar = ax.figure.colorbar(im,ax=ax)
ax.set_yticks(np.arange(len(labels)), labels=labels)
labels.reverse()
ax.set_xticks(np.arange(len(labels)+1), labels=labels+[1])
fig.savefig('heatbranch')
vel = out[1].get_velocity()

# bar.move(vel)
# bar.plot_struct('program/rigid_kin/5_bar_moved',1)
# print(regions.get_constraint_matrix())

bar.session.terminate()