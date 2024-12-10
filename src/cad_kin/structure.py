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
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
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
            constraints = self.compile_constraints()
            param_rules = []
            for elem in self.elements:
                if elem.b_parametric:
                    param_rules+=elem.param_rule
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

    def draw(self,alpha,hinge_size,params,**kwargs):
        drawing_thickness = hinge_size

        patches = []
        nodes = PatchCollection(patches,match_original=True)
        patches
        for elem in self.elements:
            p_i = params[elem.param_ids]
            patches.append(elem.plot(self.nodes,drawing_thickness,params=p_i,**kwargs))

        elem_patches = PatchCollection(patches,match_original=True)

        patches = []
        for n in self.nodes:
            c=plt.Circle((n.pos[0],n.pos[1]),drawing_thickness/2,facecolor='white',edgecolor='black',alpha=alpha,linewidth=2)
            patches.append(c)

        node_patches = PatchCollection(patches,match_original=True)

        return node_patches,elem_patches
    
    def plot(self,ax,param_vals,alpha,hinge_size,color=None,annotate=False):

        # scale = min(6.4/bar.x_dim,4/bar.y_dim)
        # print(scale)            
        
        ax.axis('off')
        # ax.set_xlim(self.x_dim*-0.2,self.x_dim*1.3)
        # ax.set_ylim(self.y_dim*-0.2,self.y_dim*1.3)
        ax.axes.set_aspect('equal')

        node, elem = self.draw(alpha,hinge_size)
        c1 = self.plot_bc(alpha)
        if isinstance(param_vals,np.ndarray):

            c2 = self.plot_params(param_vals,alpha)
        nodes,bars = self.plot_struct(alpha)
        if isinstance(color,np.ndarray):
            c1.set_facecolor(color)
            c1.set_alpha(alpha)
            if isinstance(param_vals,np.ndarray):
                c2.set_facecolor(color)
                c2.set_alpha(alpha)
            bars.set_color(color)       
            bars.set_alpha(alpha)
            nodes.set_facecolor(color)       
            # nodes.set_alpha(alpha)
        
        ax.add_collection(bars)
        if isinstance(param_vals,np.ndarray):
            c2.set_edgecolor(None)
            ax.add_collection(c2)
        
        c1.set_edgecolor(None)

        ax.add_collection(c1)
        ax.add_collection(nodes)

        # t = ax.transData.inverted()+transforms.Affine2D().scale(scale)+transforms.Affine2D().translate(1.6,4)+fig.dpi_scale_trans
        # c1.set_transform(c1.get_transform()+t)

    def plot_annotation(self,ax,param_vals):
        t = self.d0*.99

        for i, (d,p) in enumerate(zip(self.nodes_dofs,self.nodes_posns)):
            ax.text(p[0],p[1]+t,f"($v_{{{d[0]}}}$, $v_{{{d[1]}}}$)")
            ax.text(p[0]+t,p[1],f"$node_{i}$")
        for i, pval in enumerate(param_vals):
            ax.text(1,1,f"$\\alpha_{i}$",color="cyan")

    def move(self,flex):
        for i,n in enumerate(self.nodes):
            n.pos+=flex[n.dof]

        
