from cad_kin.rigidity_mech import RigidMech
from cad_kin.parametric_constraint import Parameter, ParametricConstraint

import numpy as np
import matplotlib.pyplot as plt

class Pin(RigidMech):
    eq_symbol = "=="
    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)

    def __call__(self,nodes):
        node = nodes[self.node_ids][0]

        map = node.get_map(self.n_dof)
        p = Parameter("1")
        constr = ParametricConstraint({p:map})
        return constr
    
    def get_constraint_strings(self, nodes,mod_mat):
        if self.b_parametric:
            parameter = Parameter(f"a{self.n_params}")
            # define factors for polynomial terms
            param_const = self(nodes)@mod_mat
            param_const2 = ParametricConstraint(
                {parameter:1}
            )
            param_const = param_const2*param_const

            # save information for post processing
            self.param_ids = [self.n_params]
            self.param_rule = ["bin"]

            # Incrememt Parameter Counter
            RigidMech.n_params+=1

            return super().get_constraint_strings(param_const)
        else:
            constants = self(nodes)@mod_mat
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