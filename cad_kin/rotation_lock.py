from cad_kin.rigidity_mech import RigidMech
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
class RotationLock(RigidMech):
    eq_symbol = "=="
    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)


    def __call__(self,nodes):
        nodes = nodes[self.node_ids]
        pos, dofs = self.get_node_info(nodes)

        # check current angle 
        a = pos[2:4]-pos[0:2]
        b = pos[4:6]-pos[0:2]

        cos_angle = np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b))

        # rotate angle to prevent overlap only if rotation doesn't cause overlap
        if abs(cos_angle-np.sqrt(2)/2)>1e-5:
            ref = np.array([
                [np.sqrt(2)/2,-np.sqrt(2)/2],
                [np.sqrt(2)/2,np.sqrt(2)/2],
            ])
        else:
            ref = np.identity(2)

        # change angle between elements so its not zero
        a_ref = np.matmul(ref,a)

        # precompute jacobians
        da = self.get_map_matrix(dofs[2:4])-self.get_map_matrix(dofs[0:2])
        db = self.get_map_matrix(dofs[4:6])-self.get_map_matrix(dofs[0:2])
        da_ref = np.matmul(ref,da)

        # precompute norms
        a_norm = np.linalg.norm(a)
        b_norm = np.linalg.norm(b)
        a_ref_norm = np.linalg.norm(a_ref)

        # put vectors in correct dimensions
        a_t = a[None,:]
        a = a[:,None]
        a_ref_t = a_ref[None,:]
        a_ref = a_ref[:,None]
        b_t = b[None,:]
        b = b[:,None]

        # compute constraint
        
        constr = (
            np.matmul(b_t,da_ref)+np.matmul(a_ref_t,db)
        )-(
            np.matmul(a_ref_t,b)*(np.matmul(a_ref_t,da_ref)*b_norm+np.matmul(b_t,db)*a_ref_norm)
        )/(a_ref_norm*b_norm)
        return constr
    
    def get_constraint_strings(self, nodes):

        if self.b_parametric:
            
            # define factors for polynomial terms
            param_const = self(nodes)
            
            # define map so parameters can be attributed to correct polynomial terms
            nodes = nodes[self.node_ids]
            pos, dofs = self.get_node_info(nodes)
            param_map = {
                dofs[0]:f"*a{self.n_params}",
                dofs[1]:f"*a{self.n_params}",
                dofs[2]:f"*a{self.n_params}",
                dofs[3]:f"*a{self.n_params}",
                dofs[4]:f"*a{self.n_params}",
                dofs[5]:f"*a{self.n_params}",
            }

            # save information for post processing
            self.param_ids = [self.n_params]
            self.param_rule = ["bin"]

            # Incrememt Parameter Counter
            self.n_params+=1

            return super().get_constraint_strings(param_const,[param_map])
        else:
            constants = self(nodes)
            return super().get_constraint_strings(constants)
        
    def plot(self, nodes, drawing_thickness, drawing_color='#D0D0D0',params=None):
        if self.b_parametric:
            if not params:
                return []
            if (params==0).all():
                return []
            self.angle = self.get_angle_from_params(params)

        
        pos, dofs = self.get_node_info(nodes)
        x = pos[0]
        y = pos[1]
        t = drawing_thickness*0.99
        lock = []

        a = patches.Arc(x,y,width=drawing_thickness*1.5,height=drawing_thickness*1.5,theta1=45,theta2=20)
        lock.append(a)
        a = patches.FancyArrowPatch(
            (x+np.cos(20*np.pi/180)*drawing_thickness*1.5,y+np.sin(20*np.pi/180)*drawing_thickness*1.5),
            (x+np.cos(0)*drawing_thickness*1.5,y+np.sin(0)*drawing_thickness*1.5),
            connectionstyle=f"arc3,rad={drawing_thickness*1.5}"
        )
        lock.append(a)


        if self.b_parametric:
            for s in lock:
                s.set_facecolor("#389ac7ff")

        return lock