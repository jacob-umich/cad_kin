import numpy as np

class RigidMech():
    n_params = 0
    def __init__(self,element,n_dof) -> None:
        self.node_ids = np.array(element['nodes'])
        self.n_dof = n_dof
        self.b_parametric = element.get("parametric",False)
        self.parametric_options = element.get("p_options",{})

    def get_node_info(self,nodes):
        positions = []
        for node in nodes:
            positions.append(node.pos)
        # positions = np.concatenate([node.pos for node in nodes],axis=0)
        dofs = np.concatenate([node.dof for node in nodes],axis = 0)
        return positions,dofs
    
    def get_map_matrix(self,dofs):
        map = np.zeros((2,self.n_dof))
        map[0,dofs[0]]=1
        map[1,dofs[1]]=1
        return map
    
    def get_constraint_strings(self,constants):
        precision = 10


        
        strings = constants.get_string()
        outs = []
        # each rigidity mechanisms might need multiple contraints. each might
        # have its own param map
        for s in strings:
            out = "("
            out +=s
            out+=f"){self.eq_symbol}0"
            if out == "()==0":
                out = ""
            if out == "()>=0":
                out = ""
            outs.append(out)

        return outs
    
    def plot(self, nodes, drawing_thickness, drawing_color='#D0D0D0', params=None):
        if self.b_parametric:
            if not isinstance(params,np.ndarray):
                return []
            if (abs(params)<1e-5).all():
                return []
        patches = self.plot_internal(nodes, drawing_thickness, drawing_color, params)
        for s in patches:
            if self.b_parametric:
                s.set_facecolor("#389ac7ff")
        return patches