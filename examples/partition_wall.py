from cad_kin.structure import Structure
import numpy as np

d = 5
w = 5
caps = 1.5
t= 1
l = np.sqrt(w**2+d**2)
dx = t*w/d+t*w**2/(l*d)+t*d/l
dy=t
dx2 = t*w**2/(l*d)+t*d/l
w1 = 2*(w-dy/d*w)+t
total_l = 2*caps+w+dx2+t+2*dx+w1

nodes=np.array([
    [0,0],
    [0,d],
    [caps,0],
    [caps+w,d],
    [caps+w+dx2,d],
    [caps+dx,dy],
    [caps+dx+w1,dy],
    [caps+w+dx2+t,d],
    [caps+2*dx+w1,0],
    [caps+w+dx2+t+dx,d-dy],
    [caps+w+dx2+t+dx+w1,d-dy],
    [caps+2*dx+w1+t,0],
    [caps+2*dx+w1+t+dx2,0],
    [caps+w+dx2+t+2*dx+w1,d],
    [2*caps+2*dx+w1+t+dx2+w,0],
    [2*caps+w+dx2+t+2*dx+w1,d],
])

elements = [
    {
        "type":"link",
        "nodes":[0,1]
    },
    {
        "type":"link",
        "nodes":[2,3]
    },
    {
        "type":"link",
        "nodes":[1,4],
    },
    {
        "type":"link",
        "nodes":[0,8],
    },
    {
        "type":"link",
        "nodes":[4,5],
    },
    {
        "type":"link",
        "nodes":[5,6],
    },
    {
        "type":"link",
        "nodes":[6,7],
    },
    {
        "type":"link",
        "nodes":[8,9],
    },
    {
        "type":"link",
        "nodes":[9,10],
    },
    {
        "type":"link",
        "nodes":[10,11],
    },
    {
        "type":"link",
        "nodes":[12,13],
    },
    {
        "type":"link",
        "nodes":[11,14],
    },
    {
        "type":"link",
        "nodes":[7,15],
    },
    {
        "type":"link",
        "nodes":[14,15],
    },
    {
        "type":"midspan",
        "nodes":[0,8,2]
    },
    {
        "type":"midspan",
        "nodes":[1,2,3]
    },
    {
        "type":"midspan",
        "nodes":[7,15,13]
    },
    {
        "type":"midspan",
        "nodes":[11,14,12]
    },
    {
        "type":"rotationlock",
        "nodes":[0,1,8]
    },
    {
        "type":"rotationlock",
        "nodes":[14,15,12]
    },
    {
        "type":"pin",
        "nodes":[0]
    },
    {
        "type":"roller",
        "nodes":[1],
        "angle":90
    },

]
    
part_wall = {
    "nodes":nodes,
    "elements":elements
}

part_wall = Structure(struct_dict=part_wall)

cosnt = part_wall.compile_constraints()

with open("examples/part_wal_cad.txt","w") as f:
    f.write(cosnt)
