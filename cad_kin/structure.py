import json
from cad_kin.node import Node
from cad_kin.contact_boundary import ContactBC
from cad_kin.midspan_connect import MidspanConnect
from cad_kin.rigid_link import RigidLink
from cad_kin.pin import Pin
from cad_kin.roller import Roller
from cad_kin.rotation_lock import RotationLock
from cad_kin.strut import Strut
from cad_kin.cadtree import CadTree
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wlexpr
from dotenv import load_dotenv,find_dotenv
import os
import numpy as np

class Structure():

    element_dict = {
        "link":RigidLink,
        "contactbc":ContactBC,
        "midspan":MidspanConnect,
        "pin":Pin,
        "roller":Roller,
        "strut":Strut,
        "rotationlock":RotationLock,

    }

    def __init__(self,struct_dict=None):
        if struct_dict:
            self.load_struct_dict(struct_dict)
        try:
            load_dotenv()
            self.session = WolframLanguageSession(os.getenv("WOLFRAM_KERNEL_PATH"))
        except Exception:
            print('wolfram kernel not initialized')
            
    
    def load(self,fp):
        with open(fp,"r") as f:
            struct_dict = json.load(f)
        self.load_struct_dict(struct_dict)

    def load_struct_dict(self, struct_dict):
        node_data = struct_dict["nodes"]
        self.n_dof = len(node_data)*2
        

        self.nodes = np.array([Node(data) for data in node_data])

        elem_data = struct_dict["elements"]
        self.elements = []
        for elem in elem_data:
            elem_obj = self.element_dict[elem["type"]](elem,self.n_dof)
            self.elements.append(
                elem_obj
            )
        self.n_params = RigidLink.n_params
        
    def compile_constraints(self):
        constraints = "out=CylindricalDecomposition[\n{"
        for element in self.elements:
            strings = element.get_constraint_strings(self.nodes)
            
            constraints+= ",\n".join(strings)
            if not (element==self.elements[-1]):
                constraints+=",\n"
            if isinstance(element,RigidLink) and not isinstance(element,MidspanConnect):
                strings = element.get_contact_constraint_strings(self.nodes)
                if ",\n".join(strings)=="":
                    continue
                constraints+=",\n".join(strings)
                if not (element==self.elements[-1]):
                    constraints+=",\n"
        constraints +="},\n{"

        # for k in parameters:
        #     out+=f"{k}, "
        for i in range(self.n_dof):
            if i!=self.n_dof-1:
                constraints+=f"v{self.n_dof-1-i}, "
            else:
                constraints+=f"v{self.n_dof-1-i}"+"}\n"
        constraints+=']'

        return constraints
    
    def cad(self) -> CadTree:
        try:
            param_rules = []
            for elem in self.elements:
                if elem.b_parametric:
                    param_rules+=elem.param_rule
            constraints = self.compile_constraints()
            regions =  self.session.evaluate(wlexpr(constraints))
            tree = CadTree(regions,self.n_dof,self.n_params,param_rules)
            return tree
        except Exception as e:
            print(e)
            print('wolfram kernel not initialized')       
    
    def get_labels(self): 
        dofs = []
        for i in range(self.n_dof):
            dofs.append(f"v{self.n_dof-1-i}")
        params = []
        for i in range(self.n_params):
            params.append(f"c{i}")
        return params+dofs


        


        
