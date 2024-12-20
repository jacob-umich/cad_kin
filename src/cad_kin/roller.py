from cad_kin.rigidity_mech import RigidMech
import matplotlib.pyplot as plt
import matplotlib as mpl
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
                RigidMech.n_params+=2

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
                RigidMech.n_params+=1

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
                RigidMech.n_params+=1
            
            return super().get_constraint_strings(param_const,[param_map])

        else:
            constants = self(nodes)
            return super().get_constraint_strings(constants)
        
    def get_angle_from_params(self,params):
        if self.roll_direction=="any":
            if params[1]==0:
                angle=90
            else:
                angle = np.arctan(-params[0]/params[1])*180/np.pi

        if self.roll_direction=="x":
            angle = 0

        if self.roll_direction=="y":
            angle =90

        return angle
 
    def plot_internal(self,nodes,drawing_thickness,drawing_color ='#D0D0D0',params = None):
        if self.b_parametric:
            self.angle = self.get_angle_from_params(params)

        nodes = nodes[self.node_ids]
        pos, dofs = self.get_node_info(nodes)
        x = pos[0]
        y = pos[1]
        t = drawing_thickness*0.99
        roller = []
        r = plt.Rectangle((x-0.7*t,y-t*0.7),width=t*1.4,height=t*0.7,facecolor=drawing_color)
        roller.append(r)
        c = plt.Circle((x,y),t*.7,facecolor=drawing_color)
        
        roll_1 = plt.Circle((x-0.3*t,y-t*0.9),t*.2,facecolor=drawing_color)

        roll_2 = plt.Circle((x+t*0.3,y-t*0.9),t*.2,facecolor=drawing_color) 

        base = plt.Rectangle((x-0.7*t,y-t*1.2),width=t*1.4,height=t*0.1,facecolor=drawing_color)
        roller.append(c)
        roller.append(base)
        roller.append(roll_1)
        roller.append(roll_2)

        for i in range(5):
            start = (i-2)*0.3-0.2
            mark = plt.Rectangle((x+start*t,y-t*1.5),width=t*.42,height=t*0.06,angle = 45, facecolor=drawing_color)
            roller.append(mark)

        for s in roller:
            transf = mpl.transforms.Affine2D().rotate_deg_around(x,y,self.angle)
            s.set_transform(transf)

        return roller