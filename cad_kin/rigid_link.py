from cad_kin.rigidity_mech import RigidMech
import numpy as np
from cad_kin.util import decompose

class RigidLink(RigidMech):
    eq_symbol = "=="
    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)
        self.t= element.get("thickness",1)


    def __call__(self,nodes):
        nodes = nodes[self.node_ids]

        param_nodes = [node.b_linear_parametric for node in nodes]

        # if all(param_nodes):
        #     da = nodes[0].node_j.get_map()-nodes[1].node_j.get_map()
        #     db = 

        positions, dofs = self.get_node_info(nodes)

        # vectore between two nodes
        a = positions[2:4]-positions[0:2]

        # difference in map vectors
        da = self.get_map_matrix(dofs[2:4])-self.get_map_matrix(dofs[0:2])

        constr = np.matmul(a[None,:],da)
        return constr
    
    def detect_contact(self,nodes):
        self_nodes = nodes[self.node_ids]
        contact_nodes = nodes[[not (node in self_nodes) for node in nodes]]

        elem_pos,_ = self.get_node_info(self_nodes)
        a = elem_pos[2:4]-elem_pos[:2]

        contacts = []
        for node in contact_nodes:

            b = node.pos-elem_pos[:2]

            d,gamma = decompose(a,b)
            dist = np.linalg.norm(d)
            cond1 = (dist-self.t)<1e-3
            cond2 =  np.linalg.norm(gamma)<np.linalg.norm(a) and np.dot(gamma[:,0],a)>0
            if (dist<1e-8) and cond2:
                raise Exception(
'''
Admissible direction cannot be determined automatically 
check that contacting nodes aren't in line with linkages.
'''
                )
        if cond1 and cond2:
            contacts.append(True)
        else:
            contacts.append(False)
        return contacts
    
    def detect_node_contact(self,self_node,candidate_nodes):

        contacts = []
        for node in candidate_nodes:                
            d = node.pos-self_node.pos
            cond = np.linalg.norm(d)-self.t<1e-3
            contacts.append(cond)
        return contacts
    
    def get_midspan_contact_constraint(self,self_nodes,node):
        elem_pos, dofs = self.get_node_info(self_nodes)
        a = elem_pos[2:4]-elem_pos[:2]
        b = node.pos-elem_pos[:2]
        da = self.get_map_matrix(dofs[2:4])-self.get_map_matrix(dofs[0:2])
        db = self.get_map_matrix(node.dof)-self.get_map_matrix(dofs[0:2])


        q = a/np.linalg.norm(a)
        p = q*np.dot(q,b)  
        dq = (np.linalg.norm(a)*da-np.matmul(a[:,None],np.matmul(a[None,:],da))/np.linalg.norm(a))/np.linalg.norm(a)**2
        dp = dq*np.dot(q,b)+np.matmul(q[:,None],(np.matmul(b[None,:],dq)+np.matmul(q[None,:],db)))
        constraint = np.matmul((b-p)[None,:],(db-dp))
        return constraint

    def get_node_contact_constraint(self,self_node, node):
        a = node.pos-self_node.pos
        da = self.get_map_matrix(node.dof)-self.get_map_matrix(self_node.dof)
        return np.matmul(a[None,:],da)/np.linalg.norm(a)
    
    def get_contact_constraints(self,nodes):
        constraints = []
        self_nodes = nodes[self.node_ids]

        # get all nodes not a part of the linkage
        candidate_nodes = nodes[[not (node in self_nodes) for node in nodes]]
        
        # select nodes in contact with linkage midspan
        mask = self.detect_contact(nodes)
        contact_nodes = candidate_nodes[mask]
        for node in contact_nodes:
            constraints.append(self.get_midspan_contact_constraint(self_nodes,node))

        # Select nodes that could still be in contact with linkage nodes
        candidate_nodes = candidate_nodes[not(mask)]

        # generate contact constraints for ith node of linkage
        mask1 = self.detect_node_contact(self_nodes[0],candidate_nodes)
        contact_nodes = candidate_nodes[mask]
        for node in contact_nodes:
            constraints.append(self.get_node_contact_constraint(self_nodes[0],node))

        # generate contact constraints for jth node of linkage
        mask2 = self.detect_node_contact(self_nodes[1],candidate_nodes)
        contact_nodes = candidate_nodes[mask]
        for node in contact_nodes:
            constraints.append(self.get_node_contact_constraint(self_nodes[1],node))

        return constraints


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