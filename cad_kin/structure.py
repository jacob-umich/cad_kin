import json
from cad_kin.node import Node
from cad_kin.contact_boundary import ContactBC
from cad_kin.midspan_connect import MidspanConnect
from cad_kin.rigid_link import RigidLink
from cad_kin.pin import Pin
from cad_kin.roller import Roller
from cad_kin.rotation_lock import RotationLock
from cad_kin.strut import Strut
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wlexpr
from dotenv import load_dotenv

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
            self.session = WolframLanguageSession(os.getenv("WOLFRAM_KERNEL_PATH"),kernel_loglevel=logging.DEBUG)
        except Exception:
            print('wolfram kernel not initialized')
            
    
    def load(self,fp):
        with open(fp,"r") as f:
            struct_dict = json.load(f)
        self.load_struct_dict(struct_dict)

    def load_struct_dict(self, struct_dict):
        node_data = struct_dict["nodes"]
        self.n_dof = len(node_data)*2

        self.nodes = [Node(data) for data in node_data]

        elem_data = struct_dict["elements"]
        self.elements = []
        for elem in elem_data:
            self.elements.append(
                self.element_dict[elem["type"]](elem,self.n_dof)
            )
        
    def compile_constraints(self):
        constraints = "out=CylindricalDecomposition[\n{"
        for element in self.elements:
            constraints+=element.get_constraint_string(self.nodes)
            if not (element==self.elements[-1]):
                constraints+=",\n"
            if isinstance(element,RigidLink):
                constraints+=element.get_contact_constraint_strings(self.nodes)
                if not (element==self.elements[-1]):
                    constraints+=",\n"
        constraints +="},\n{"

        # for k in parameters:
        #     out+=f"{k}, "
        for i in range(self.n_dof):
            if i!=self.n_dof-1:
                constraints+=f"x{self.n_dof-1-i}, "
            else:
                constraints+=f"x{self.n_dof-1-i}"+"}\n"
        constraints+=']'
        self.constraints = constraints
        return constraints
    
    def cad(self):
        try:
            return self.session.evaluate(wlexpr(self.constraints))
        except Exception as e:
            print(e)
            print('wolfram kernel not initialized')


        


        
