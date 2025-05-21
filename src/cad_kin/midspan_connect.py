from cad_kin.rigid_link import RigidLink
from cad_kin.rotation_lock import RotationLock
from cad_kin.rigidity_mech import RigidMech
import numpy as np
from copy import copy

class MidspanConnect(RigidLink,RotationLock):
    eq_symbol = "=="

    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)

    def __call__(self,nodes):
        const1 = RigidLink.__call__(self,nodes)
        const2 = RotationLock.__call__(self,nodes)

        return const1.concatenate(const2,axis=0 )
    
    def get_constraint_strings(self, nodes,mod_mat):

        if self.b_parametric:
            return ""
        else:
            constants = self(nodes)@mod_mat
            return super(RotationLock,self).get_constraint_strings(constants)
        
