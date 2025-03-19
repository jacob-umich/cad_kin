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

        return np.concatenate([const1, const2],axis=0)
    
    def get_constraint_strings(self, nodes,mod_mat):

        if self.b_parametric:
            # i am pretty sure constraint is no longer needed if constraint is parametric
            # replace node for parametric definition

            # param_node = nodes[self.node_ids[1]]
            
            # param_node = param_node.make_linear_parametric(
            #     nodes[self.node_ids[0]],
            #     nodes[self.node_ids[2]],
            #     self.n_params
            # )

            # nodes[self.node_ids[1]]=param_node

            # # define factors for polynomial terms
            # param_const = np.matmul(self(nodes),mod_mat)
            
            # # define map so parameters can be attributed to correct polynomial terms
            # nodes = nodes[self.node_ids]
            # pos, dofs = self.get_node_info(nodes)
            # polynomial1 = f"*(1-2*a{self.n_params}+a{self.n_params}^2)"
            # param_map1 = {
            #     dofs[0]:polynomial1,
            #     dofs[1]:polynomial1,
            #     dofs[4]:polynomial1,
            #     dofs[5]:polynomial1,
            # }
            # polynomial2 = f"*(1-3*a{self.n_params}+3*a{self.n_params}^2-a{self.n_params}^3)"
            # param_map2 = {
            #     dofs[0]:polynomial2,
            #     dofs[1]:polynomial2,
            #     dofs[4]:polynomial2,
            #     dofs[5]:polynomial2,
            # }
            # param_maps = [param_map1,param_map2]

            # # save information for post processing
            # self.param_ids = [self.n_params]
            # self.param_rule = ["bin"]

            # # Incrememt Parameter Counter
            # # self.n_params+=1

            # # add custom constraints to bound parameter
            # result = (
            #     super(RotationLock,self).get_constraint_strings(param_const,param_maps)
            #     +[
            #         f"a{self.n_params}>=0",
            #         f"a{self.n_params}<=1"
            #     ]
            # )

            return ""
        else:
            constants = np.matmul(self(nodes),mod_mat)
            return super(RotationLock,self).get_constraint_strings(constants)
        
