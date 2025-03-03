from cad_kin.rigidity_mech import RigidMech
import numpy as np
import matplotlib.pyplot as plt

class Pin(RigidMech):
    eq_symbol = "=="
    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)

    def __call__(self,nodes):
        node = nodes[self.node_ids][0]

        map = self.get_map_matrix(node.dof)
        return map
    
    def get_constraint_strings(self, nodes,mod_mat):
        if self.b_parametric:
            
            # define factors for polynomial terms
            param_const = np.matmul(self(nodes),mod_mat)
            
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
            RigidMech.n_params+=1

            return super().get_constraint_strings(param_const,[param_map]*2)
        else:
            constants = np.matmul(self(nodes),mod_mat)
            return super().get_constraint_strings(constants)
    
    def plot(self,nodes,drawing_thickness,drawing_color ='#D0D0D0',params = None ):

        pos, dofs = self.get_node_info(nodes)
        x = pos[0]
        y = pos[1]
        t = drawing_thickness*0.99
        hinge = []
        r = plt.Rectangle((x-0.7*t,y-t*0.6),width=t*1.4,height=t*0.6,facecolor=drawing_color)
        hinge.append(r)
        c = plt.Circle((x,y),t*.7,facecolor=drawing_color) 
        base = plt.Rectangle((x-0.9*t,y-t*0.7),width=t*1.8,height=t*0.1,facecolor=drawing_color)
        hinge.append(c)
        hinge.append(base)

        for i in range(5):
            start = (i-2)*0.3-0.2
            mark = plt.Rectangle((x+start*t,y-t),width=t*.42,height=t*0.06,angle = 45, facecolor=drawing_color)
            hinge.append(mark)

        return hinge