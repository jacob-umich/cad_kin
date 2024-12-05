# Cylindrical Algebraic Decomposition of Linkage Kinematics

linkages are an important idealization for many mechanisms. This library provides the tools to formulate the kinematic constraints of a linkage system and apply the CAD analysis to them.

# Requirements

# Creating a linkage to analyze:
To create a linkage, first a structure object must be created. If a dictionary has been defined in the format of the [input file JSON](#general-input-file-structure), it can be used to instantiate the structure:

    from cad_kin.structure import Structure
    s = Structure(structure_dict=None)

If a structure dict was not provided and you are creating a structure from an [input file](#general-input-file-structure), the input file can be loaded

    s.load(file_path)

The structure is now loaded. To perform the CAD algorithm and decompose the design space, use the `cad()` method, which returns a `CADTree` object 

    tree = s.cad()



# General Input File Structure:
The structure is made up of two items: nodes to define the geometry, and "elements" to define constraints between nodes. Each of these items are lists.

```
{
    "nodes":[
        ...
    ],
    "elements":[
        ...
    ],
    
}
```


## Nodes:
defines the end points of linkages and acts as a connection point between linkages and other 1-dimensional elements

json structure:
```
{
    "nodes":[
        [x,y], \\ position of node in 2-D
        ...
    ],
    ...
}
```

Node ids are assigned incrementally based on their order in the input file starting with 0.

data path: data["nodes"]
data type: list of float lists representing coordinate of each node.

## Elements and Boundary Conditions

The elements are contained in the main element block
Types:
- Rigid link Element
- Strut/Cable Element
- Midspan Linkage Connection
- Rotational Hinge Lock
- Roller Boundary Condition
- Pinned Boundary Condition
- Boundary Contact Condition

```
{
    "nodes":[
        ...
    ],
    "elements":[
        {
            \\element key-value pairs
        },
        ...
    ]
}
```
data path: data["elements"]
data type: list of objects that define constraint information.

### Rigid link element
An element rigidly connecting two nodes. This element is designated by `"type":"link"` key-value
pair 
```
{
    "nodes":[
        ...
    ],
    "elements":[
        {
            "type":"link",
            "nodes":[Number,Number]
            "thickness": Number,
            "parametric": Bool,
        },
        ...
    ]
}
```

Key-value pairs:
- nodes: The node number of beginning and end nodes. Node id's are assigned incrementally based on their order in the input file
- thickness: (Optional) Designates the thickness of the link for self-contact purposes. The width of the link is centered on its 1D definition. The default value is `1`
- parametric: (Optional) Determines whether this element will vary in the design space. The default value is `False`

### Midspan Connect
Specifies that a given node is rigidly connected to the mispan defined by two separate nodes. Nodes must all be colinear.  This element is designated by `"type":"midspan"` key-value pair. This type of constriant does not currently have the option to be parameterized in the design.

```
{
    "nodes":[
        ...
    ],
    "elements":[
        {
            "type":"midspan",
            "nodes":[Number,Number,Number]
        },
        ...
    ]
}
```

Key-value pairs:
- nodes: The node number of beginning, midspan, and end nodes, in that order. Node id's are assigned incrementally based on their order in the input file.


### Roller
Creates a roller boundary condition at the specified node. This element is designated by `"type":"roller"` key-value pair. 

```
{
    "nodes":[
        ...
    ],
    "elements":[
        {
            "type":"roller",
            "nodes":[Number],
            "angle": Number,
            "parametric": Bool,
            "p_options": {
                "roll_direction": str["any"|"x"|"y"]
            }
        },
        ...
    ]
}
```

Key-value pairs:
- nodes: The node number where the boundary will be applied. Node id's are assigned incrementally based on their order in the input file.
- angle: (Optional) The angle of the roller from horizontal in degrees. For example, if the angle is 0, the node can freely translate in the x-direction. The default value is `0`
- parametric: (Optional) Determines whether this element will vary in the design space. The default value is `False`
- p_options: (Optional) A dictionary that determines additional information needed to parameterize the roller. Key-value pair for a roller is `"roll_direction"` paired with `"any","x","y"`. If set to any, the design space will include the roller at any direction. If set to `"x"` or `"y"`, the design space will only consider the roller allowing x or y translation, respectively.

### Pin
Creates a pinned boundary condition at the specified node. This element is designated by `"type":"pin"` key-value pair. 

```
{
    "nodes":[
        ...
    ],
    "elements":[
        {
            "type":"pin",
            "nodes":[Number],
            "parametric": Bool,
        },
        ...
    ]
}
```

Key-value pairs:
- nodes: The node number where the boundary will be applied. Node id's are assigned incrementally based on their order in the input file.
- parametric: (Optional) Determines whether this element will vary in the design space. The default value is `False`.

### Strut/cable element
An element connecting two nodes. Struts prevent the connected nodes moving closer, and cables prevent them moving apart. These elements are designated by `"type":"strut"` key-value
pair 
```
{
    "nodes":[
        ...
    ],
    "elements":[
        {
            "type":"strut",
            "nodes":[Number,Number],
            "cable":Bool,
            "thickness": Number,
            "parametric": Bool,
        },
        ...
    ]
}
```

Key-value pairs:
- nodes: The node number of beginning and end nodes. Node id's are assigned incrementally based on their order in the input file
- cable: (Optional) When set to `True`, converts strut into a cable. The default value is `False`.
- thickness: (Optional) Designates the thickness of the link for self-contact purposes. The width of the link is centered on its 1D definition. The default value is `1`.
- parametric: (Optional) Determines whether this element will vary in the design space. The default value is `False`.

### Rotation Lock
Locks the angle between two connected links. This element is designated by `"type":"rotationlock"` key-value pair 

```
{
    "nodes":[
        ...
    ],
    "elements":[
        {
            "type":"rotationlock",
            "nodes":[Number,Number,Number],
            "parametric": Bool,
        },
        ...
    ]
}
```

Key-value pairs:
- nodes: The node number of the two links being locked. The first is from one link, the second is the shared node, and the third is from the other link. Node id's are assigned incrementally based on their order in the input file.

- parametric: (Optional) Determines whether this element will vary in the design space. The default value is `False`.

### Boundary Contact Condition
Creates a boundary contact condition at the specified node. This element is designated by `"type":"contactbc"` key-value pair. 

```
{
    "nodes":[
        ...
    ],
    "elements":[
        {
            "type":"contactbc",
            "nodes":[Number],
            "angle": Number,
            "parametric": Bool,
            "p_options": {
                "roll_direction": str["any"|"x"|"y"]
            }
        },
        ...
    ]
}
```

Key-value pairs:
- nodes: The node number where the boundary will be applied. Node id's are assigned incrementally based on their order in the input file.
- angle: (Optional) The angle of the contact surface from horizontal in degrees. For example, if the `"angle": 0`, the node can freely translate in the x-direction and in the postive vertical direction. The default value is `0`
- parametric: (Optional) Determines whether this element will vary in the design space. The default value is `False`
- p_options: (Optional) A dictionary that determines additional information needed to parameterize the boundary. Key-value pair for a boundary contact is `"roll_direction"` paired with `"any","x","y"`. If set to any, the design space will consider the contact surface at any direction. If set to `"x"` or `"y"`, the design space will only consider the surface to be horizontal or vertical, respectively.