import numpy as np

class RigidMech():
    n_params = 0
    def __init__(self,element,n_dof) -> None:
        self.node_ids = element['nodes']
        self.n_dof = n_dof
        self.b_parametric = element.get("parametric",False)
        self.parametric_options = element.get("p_options",{})

    def get_node_info(self,nodes):
        positions = np.concatenate([node.pos for node in nodes],axis=0)
        dofs = np.concatenate([node.dof for node in nodes],axis = 0)
        return positions,dofs
    
    def get_map_matrix(self,dofs):
        map = np.zeros((2,self.n_dof))
        map[0,dofs[0]]=1
        map[1,dofs[1]]=1
        return map
    
    def get_constraint_strings(self,constants,param_maps=None):
        precision = 5

        if not param_maps:
            param_maps = [{}]*len(constants)
        
        strings = []
        # each rigidity mechanisms might need multiple contraints. each might
        # have its own param map
        for constant,param_map in zip(constants,param_maps):
            if not any(constant):
                continue
            out = "("
            first=True
            for i,v in enumerate(constant):
                if abs(v)>1e-15:
                    param_string = f"{param_map.get(i,'')}"
                    if v<0:
                        s = "-"
                    else:
                        s = "+"
                    if( first and s=="+"):
                        out+=f" {int(abs(v)*10**(precision))}*v{i}{param_string} "
                        first=False
                    else:
                        out+=f" {s}{int(abs(v)*10**(precision))}*v{i}{param_string} "
                        first = False
            out+=f"){self.eq_symbol}0"
            strings.append(out)

        return strings
    
    def plot(self, nodes, drawing_thickness, drawing_color='#D0D0D0', params=None):
        return []