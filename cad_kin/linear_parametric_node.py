from cad_kin.node import Node

class LinearParametricNode(Node):
    b_linear_paramtric = True
    def __init__(self,original_node,node_i,node_j,parameter):
        self.id = original_node.id
        self.dof = original_node.dof
        self.node_i = node_i
        self.node_j = node_j
        self.parameter = parameter

    def get_map(self):
        b_map = self.node_j.get_map()
        diff_map= self.node_i.get_map()-self.node_j.get_map()
        return b_map, diff_map

    