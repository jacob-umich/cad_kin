from cad_kin.roller import Roller
import numpy as np

class ContactBC(Roller):
    eq_symbol = ">="
    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)

    def __call__(self,nodes):
        return super().__call__(nodes)

    