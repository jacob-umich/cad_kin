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
from cad_kin.rigidity_mech import RigidMech
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wlexpr
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from dotenv import load_dotenv,find_dotenv
import os
import numpy as np
import scipy.linalg
import traceback
import logging

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
        load_dotenv()
        self.wolfram_path = os.getenv("WOLFRAM_KERNEL_PATH")

            
    
    def load(self,fp):
        with open(fp,"r") as f:
            struct_dict = json.load(f)
        self.load_struct_dict(struct_dict)

    def load_struct_dict(self, struct_dict):
        node_data = struct_dict["nodes"]
        self.n_dof = len(node_data)*2
        
        # will need to change how DOF are defined if we have parametric nodes
        self.nodes = np.array([Node(data) for data in node_data])

        # find bounds of structure
        x_min = min([node.pos[0] for node in self.nodes])
        x_max = max([node.pos[0] for node in self.nodes])
        y_min = min([node.pos[1] for node in self.nodes])
        y_max = max([node.pos[1] for node in self.nodes])
        self.bounds = np.array([x_min,x_max,y_min,y_max])


        elem_data = struct_dict["elements"]
        self.elements = []
        for elem in elem_data:
            elem_obj = self.element_dict[elem["type"]](elem,self.n_dof)
            self.elements.append(
                elem_obj
            )
            # set parametric node
        
    def get_modes(self):
        constraint_matrix = []
        for element in self.elements:
            if element.eq_symbol=="==" and (not element.b_parametric):
                constraint = element(self.nodes)
                constraint_matrix.append(constraint)

        constraint_matrix = np.concatenate(constraint_matrix,axis=0)
        rank = np.linalg.matrix_rank(constraint_matrix)
        q,r,p=scipy.linalg.qr(
            np.transpose(constraint_matrix),
            mode="full",
            pivoting=True
        )
        self.modes = q[:,rank:]
        return q[:,rank:]

    def compile_constraints(self,b_spectral=False):
        # bug, if compile constraints is called more than once, 
        # design parameter numbers will be incremented again
        constraints = "out=CylindricalDecomposition[\n{"

        if b_spectral:
            mod_mat = self.get_modes()
            n_dof = len(mod_mat[0])
        else:
            mod_mat = np.identity(self.n_dof)
            n_dof = self.n_dof
        for element in self.elements:
            strings = element.get_constraint_strings(self.nodes,mod_mat)
            
            
            constraints+= ",\n".join(strings)
            if not (element==self.elements[-1]):
                constraints+=",\n"
                
            if isinstance(element,RigidLink) and not isinstance(element,MidspanConnect):
                strings = element.get_contact_constraint_strings(self.nodes,mod_mat)
                if ",\n".join(strings)=="":
                    continue
                constraints+=",\n".join(strings)
                if not (element==self.elements[-1]):
                    constraints+=",\n"
        constraints +="},\n{"
        constraints = constraints.replace("0==0,\n","")

        self.n_params = RigidMech.n_params
        
        for k in range(self.n_params):
            constraints+=f"a{k}, "
        for i in range(n_dof):
            if i!=n_dof-1:
                constraints+=f"v{n_dof-1-i}, "
            else:
                constraints+=f"v{n_dof-1-i}"+"}\n"
        constraints+=']'

        return constraints
    
    def cad(self,b_spectral=False,b_debug=False) -> CadTree:
        try:
            if b_debug:
                self.session = WolframLanguageSession(self.wolfram_path,kernel_loglevel=logging.DEBUG)
                logging.basicConfig(level=logging.DEBUG)
            else:
                self.session = WolframLanguageSession(self.wolfram_path)

        except Exception as e:
            print('wolfram kernel not initialized') 
            print(e)
            return None
        try:
            constraints = self.compile_constraints(b_spectral)
            logging.log(logging.DEBUG,constraints)
            param_rules = []
            for elem in self.elements:
                if elem.b_parametric:
                    param_rules+=elem.param_rule
            regions =  self.session.evaluate(wlexpr(constraints))

            tree = CadTree(regions,self.n_dof,self.n_params,self.get_labels(),param_rules)
            return tree
        except Exception as e:
            print("CAD algorithm Failed")
            print(e)
            print(traceback.format_exc())
        finally:
            self.session.terminate()
    
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

        elem_patches = []
        for elem in self.elements:
            if elem.b_parametric:
                print(elem.param_ids)
                p_i = params[elem.param_ids]
                elem_patches.append(
                    PatchCollection(
                        elem.plot(self.nodes,drawing_thickness,params=p_i,**kwargs),
                        match_original=True
                    )
                )
            else:
                elem_patches.append(
                    PatchCollection(
                        elem.plot(self.nodes,drawing_thickness,**kwargs),
                        match_original=True
                    )
                )

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
        ax.set_xlim(
            self.bounds[0]-hinge_size*3,
            self.bounds[1]+hinge_size*3
        )
        ax.set_ylim(
            self.bounds[2]-hinge_size*3,
            self.bounds[3]+hinge_size*3
        )
        ax.axes.set_aspect('equal')

        node, elem = self.draw(alpha,hinge_size,param_vals)


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
        
        for elem_patch in elem:

            ax.add_collection(elem_patch)

        ax.add_collection(node)
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

        
