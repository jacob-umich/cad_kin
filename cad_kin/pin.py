from cad_kin.rigidity_mech import RigidMech
import numpy as np

class Pin(RigidMech):
    eq_symbol = "=="
    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)

    def __call__(self,nodes):
        node = nodes[self.node_ids][0]

        map = self.get_map_matrix(node.dof)
        return map
    
    def get_constraint_strings(self, nodes):
        if self.b_parametric:
            
            # define factors for polynomial terms
            param_const = self(nodes)
            
            # define map so parameters can be attributed to correct polynomial terms
            nodes = nodes[self.node_ids]
            pos, dofs = self.get_node_info(nodes)

            # both constraints can share the same param map
            param_map = {
                dofs[0]:f"*a{self.n_params}",
                dofs[1]:f"*a{self.n_params}",
            }

            # save information for post processing
            self.param_ids = [self.n_params]
            self.param_rule = ["bin"]

            # Incrememt Parameter Counter
            self.n_params+=1

            return super().get_constraint_strings(param_const,[param_map]*2)
        else:
            constants = self(nodes)
            return super().get_constraint_strings(constants)