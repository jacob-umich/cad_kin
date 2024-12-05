from cad_kin.rigid_link import RigidLink
import numpy as np

class Strut(RigidLink):
    eq_symbol=">="
    def __init__(self,element,n_dof):
        super().__init__(element,n_dof)
        self.b_cable = element.get("cable",False)

    def __call__(self,nodes):
        if self.b_cable:
            return -super().__call__(nodes)
        else:
            return super().__call__(nodes)
    def get_constraint_strings(self, nodes):
        string = super().get_constraint_strings(nodes)
        self.param_rule = ["flip"]
        return string
    
    def plot(self, nodes, drawing_thickness, drawing_color='#D0D0D0',params=None):
        if self.b_parametric:
            if not params:
                return []
            if (params==0).all():
                return []
            if params[0]<0:
                self.b_cable
        r = super().plot(nodes, drawing_thickness, drawing_color)
        r.set_facecolor("none")
        if self.b_cable:
            r.set_height(r.get_height/10)