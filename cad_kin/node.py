import numpy as np
from cad_kin.linear_parametric_node import LinearParametricNode

class Node():
    current_id = 0
    def __init__(self,pos):
        self.pos = np.array(pos)
        self.id = Node.current_id
        Node.current_id +=1
        self.dof = np.array([self.id*2,self.id*2+1])

    def get_map(self,n_dofs):
        map = np.zeros((2,n_dofs))
        map[0,self.dof[0]]=1
        map[1,self.dof[1]]=1
        return map
    
    def make_linear_parametric(self,node_i,node_j,parameter):
        return LinearParametricNode(self,node_i,node_j,parameter)


