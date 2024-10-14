from cad_kin.rigidity_mech import RigidMech
import numpy as np

class Roller(RigidMech):
    eq_symbol = "=="
    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)
        self.angle = element.get("angle",0)

    def __call__(self,nodes):

        node = nodes[self.node_ids][0]

        a = -np.sin(self.angle*np.pi/180)
        b = np.cos(self.angle*np.pi/180)
        values = np.array([[a,b]])
        map = self.get_map_matrix(node.dof)
        return np.matmul(values,map)
    
    def get_constraint_strings(self,nodes):
        node = nodes[self.node_ids][0]

        if self.b_parametric:
            
            roll_direction = self.parametric_options.get("roll_direction","any")
            # save roll_dir
            self.roll_direction = roll_direction

            if roll_direction=="any":

                # define map so parameters can be attributed to correct polynomial terms
                param_map = {
                    node.dof[0]:f"*a{self.n_params}",
                    node.dof[1]:f"*a{self.n_params+1}",
                }

                # define factors for polynomial terms
                self.angle = -45
                param_const = self(nodes)

                # save information for post processing
                self.param_ids = [self.n_params,self.n_params+1]
                self.param_rule = ["a","b"]

                # Incrememt Parameter Counter
                self.n_params+=2

            if roll_direction=="x":

                # define map so parameters can be attributed to correct polynomial terms
                param_map = {
                    node.dof[1]:f"*a{self.n_params}",
                }
                
                # define factors for polynomial terms
                self.angle = 0
                param_const = self(nodes)

                # save information for post processing
                self.param_ids = [self.n_params]
                self.param_rule = ["bin"]

                # Incrememt Parameter Counter
                self.n_params+=1

            if roll_direction=="y":
                # define map so parameters can be attributed to correct polynomial terms
                param_map = {
                    node.dof[0]:f"*a{self.n_params}",
                }
                
                # define factors for polynomial terms
                self.angle = -90
                param_const = self(nodes)

                # save information for post processing
                self.param_ids = [self.n_params]
                self.param_rule = ["bin"]

                # Incrememt Parameter Counter
                self.n_params+=1
            
            return super().get_constraint_strings(param_const,[param_map])

        else:
            constants = self(nodes)
            return super().get_constraint_strings(constants)