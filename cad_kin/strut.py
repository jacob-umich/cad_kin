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
