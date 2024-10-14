class ParametricConstraint():
    def __init__(self, param_dict):
        self.param_dict = param_dict

    def __matmul__ (self, other):
        b_pcosntraint = isinstance(other,ParametricConstraint)

        if b_pcosntraint:
            for param_i,const_i in self.param_dict.items():
                for param_j,const_j in other.param_dict.items():
                    # new_param = 