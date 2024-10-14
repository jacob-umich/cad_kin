# Cylindrical Algebraic Decomposition of Linkage Kinematics

linkages are an important idealization for many mechanisms. This library provides the tools to formulate the kinematic constraints of a linkage system and apply the CAD analysis to them.

# Requirements


## Nodes:
defines the end points of linkages and acts as a connection point between linkages and other 1-dimensional elements

json structure:
```
{
    "nodes":[
        [coord1,coord2] \\ position of node in 2-D
    ]
}
```

Node ids are assigned incrementally based on their order in the input file starting with 0.
data path: data["nodes"]
data type: list of float lists representing coordinate of each node

## Elements
Element types:
- link
- midspan_connect
- rotation_lock

### Rigid link element

properties:
- nodes: node number of beginning and end nodes

### Midspan Connect

3rd node defined is midspan node

2nd node needs to be midspan node to be compatible with rigid link

### Roller

specify angle in degrees

### Parametric constraints
 Parameters are automatically incremented like nodes