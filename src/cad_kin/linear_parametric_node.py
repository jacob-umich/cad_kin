from cad_kin.node import Node
from cad_kin.parametric_constraint import Parameter, ParametricConstraint
import numpy as np

class LinearParametricNode(Node):
    b_linear_paramtric = True
    def __init__(self,node_i,node_j,parameter):
        self.dof = np.concatenate([node_i.dof,node_j.dof],axis=0)
        self.node_i = node_i
        self.node_j = node_j
        self.parameter = parameter
        self.param_obj = Parameter(parameter)
        self.constant_obj = Parameter("1")
        self.pos = ParametricConstraint({
            self.param_obj:(node_j.pos-node_i.pos),
            self.constant_obj: node_i.pos
        })
    def get_map(self,n_dofs):
        b_map = self.node_i.get_map(n_dofs)
        diff_map= self.node_j.get_map(n_dofs)-self.node_i.get_map(n_dofs)
        return ParametricConstraint({
            self.param_obj:diff_map,
            self.constant_obj: b_map
        })

class ParametricNode(Node):
    b_parametric = True
    def __init__(self,node_i,parameter):
        self.dof = node_i.dof
        self.x_p = Parameter(parameter[0])
        self.y_p = Parameter(parameter[1])
        self.constant_obj = Parameter("1")
        self.pos = ParametricConstraint({
            self.x_p:np.array([[1,0]]),
            self.y_p:np.array([[0,1]]),
            self.constant_obj: node_i.pos
        })
    def get_map(self,n_dofs):
        b_map = self.node_i.get_map(n_dofs)
        diff_map= self.node_j.get_map(n_dofs)-self.node_i.get_map(n_dofs)
        return ParametricConstraint({
            self.param_obj:diff_map,
            self.constant_obj: b_map
        })