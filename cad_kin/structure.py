class Structure():
    def __init__(self,struct_dict):
        node_data = struct_dict["nodes"]
        n_dof = len(node_data)*2