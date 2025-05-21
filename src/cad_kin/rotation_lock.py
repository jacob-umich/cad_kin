from cad_kin.rigidity_mech import RigidMech
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from cad_kin.parametric_constraint import Parameter, ParametricConstraint
class RotationLock(RigidMech):
    eq_symbol = "=="
    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)


    def __call__(self,nodes)->ParametricConstraint:
        nodes = nodes[self.node_ids]
        pos, dofs = self.get_node_info(nodes)

        # check current angle 
        a = ParametricConstraint({Parameter("1"):pos[1]-pos[0]})
        b = ParametricConstraint({Parameter("1"):pos[2]-pos[0]})

        cos_angle = (a.dot(b)/np.sqrt(a.dot(a)*b.dot(b))).param_dict["1"]

        # rotate angle to prevent overlap only if rotation doesn't cause overlap
        if abs(cos_angle-np.sqrt(2)/2)>1e-5:
            ref = np.array([
                [np.sqrt(2)/2,-np.sqrt(2)/2],
                [np.sqrt(2)/2,np.sqrt(2)/2],
            ])
        else:
            ref = np.identity(2)

        # change angle between elements so its not zero
        a_ref = ref@a

        # precompute jacobians
        da = nodes[1].get_map(self.n_dof)-nodes[0].get_map(self.n_dof)
        db = nodes[2].get_map(self.n_dof)-nodes[0].get_map(self.n_dof)
        da_ref = ref@da

        # precompute norms
        a_norm = a.dot(a).sqrt()
        b_norm = b.dot(b).sqrt()
        a_ref_norm = a_ref.dot(a_ref).sqrt()

        # put vectors in correct dimensions
        a_t = a[None,:]
        a = a.transpose()
        a_ref_t = a_ref
        a_ref = a_ref.transpose()
        b_t = b
        b = b.transpose()

        # compute constraint
        
        constr = (
            b_t@da_ref+a_ref_t@db
        )-(
            a_ref_t@b*(a_ref_t@da_ref*b_norm+b_t@db*a_ref_norm)
        )/(a_ref_norm*b_norm)
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
        
    def plot(self, nodes, drawing_thickness, drawing_color='#D0D0D0',params=None):

        
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


        return lock