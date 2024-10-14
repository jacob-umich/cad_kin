import rigid_kin
import numpy as np
from contact_test import probe_contact
wall = rigid_kin.Structure('program/rigid_kin/part_wall.json')
assert(all(wall.nodes_posns[0]==[0,0]))
assert(np.linalg.norm(wall.nodes_posns[11]-np.array([15.121,0]))<0.1)

assert(wall.nodes_dofs[3,0]==6)
assert(wall.nodes_dofs[3,1]==7)


assert(wall.links[0]['pos'][3]==wall.nodes_posns[1,1])
assert(wall.straights[0]['pos'][3]==wall.nodes_posns[3,1])
assert(wall.straights[0]['pos'][11]==wall.nodes_posns[4,1])

assert(wall.n_dof==32)

wall.gen_constraints()
print(wall.constraints.shape)
print(wall.contact_constraints.shape)
# wall.plot_struct('program/rigid_kin/part_wall',1)

# bar_4 = rigid_kin.Structure('program/rigid_kin/4_bar.json')
# bar_4.gen_constraints()
# print(bar_4.constraints.shape)
# print(bar_4.contact_constraints.shape)
# bar_4.plot_struct('program/rigid_kin/test',1)
# regions = bar_4.cad('program/rigid_kin/4_bar.wls')
# bar_4.session.terminate()

# testing out constraints

#testing node to node contact
test_vector = np.zeros((wall.n_dof,1))
test_vector[8,0]=1
print(wall.contact_constraints[1,:])
dofs = [8,9]
const_n = 1
cases=[
    {
        "disp":[-1,0],
        "valid":1
    },
    {
        "disp":[1,0],
        "valid":-1
    },
    {
        "disp":[0,-1],
        "valid":1
    },
    {
        "disp":[0,1],
        "valid":1
    },
    {
        "disp":[-1,-1],
        "valid":1
    },
    {
        "disp":[1,1],
        "valid":-1
    },
    {
        "disp":[1,-1],
        "valid":-1
    },
    {
        "disp":[-1,1],
        "valid":1
    },
]
probe_contact(dofs,const_n,wall.contact_constraints,cases)

const = wall.gen_contact_constraint(wall.nodes[7],wall.links[1])
dofs = [14,15]
cases=[
    {
        "disp":[-1,0],
        "valid":-1
    },
    {
        "disp":[1,0],
        "valid":1
    },
    {
        "disp":[0,-1],
        "valid":1
    },
    {
        "disp":[0,1],
        "valid":1
    },
    {
        "disp":[-1,-1],
        "valid":-1
    },
    {
        "disp":[1,1],
        "valid":1
    },
    {
        "disp":[1,-1],
        "valid":1
    },
    {
        "disp":[-1,1],
        "valid":-1
    },
]
probe_contact(dofs,0,const,cases)
const = wall.gen_contact_constraint(wall.nodes[7],wall.links[3])
dofs = [14,15]
cases=[
    {
        "disp":[-1,0],
        "valid":-1
    },
    {
        "disp":[1,0],
        "valid":1
    },
    {
        "disp":[0,-1],
        "valid":1
    },
    {
        "disp":[0,1],
        "valid":1
    },
    {
        "disp":[-1,-1],
        "valid":-1
    },
    {
        "disp":[1,1],
        "valid":1
    },
    {
        "disp":[1,-1],
        "valid":1
    },
    {
        "disp":[-1,1],
        "valid":-1
    },
]
probe_contact(dofs,0,const,cases)