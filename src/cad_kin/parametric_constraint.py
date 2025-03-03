import numpy as np
import copy
class Parameter():
    def __init__(self, *symbols):
        if symbols[0]=="1":
            self.b_const = True
            self.ids = [-1]
            self.symbol="1"
        else:
            self.b_const=False
            self.ids = sorted([int(symbol[1:]) for symbol in symbols])
            self.symbol="*".join([f"a{id}" for id in self.ids])


    def __hash__(self):
        return self.symbol
    def __eq__(self, value):
        return self.ids == value.ids

    def __mul__(self,value):
        if self.b_const:
            # the case where this object is constant
            return value
        else:
            # the case where the other object is constant
            if value.b_const:
                return self
            else:
            # both are non constant
                return Parameter(*([f"a{id}" for id in self.ids]+[f"a{id}" for id in value.ids]))


class ParametricConstraint():
    def __init__(self, param_dict):
        self.param_dict = param_dict

    def __matmul__ (self, other):
        b_pcosntraint = isinstance(other,ParametricConstraint)

        if b_pcosntraint:
            out_dict = {}
            for param_i,const_i in self.param_dict.items():
                for param_j,const_j in other.param_dict.items():
                    key = param_i*param_j
                    values = out_dict.get(key,None)
                    new_param = np.matmul(const_i,const_j)
                    if values:
                        out_dict[key]+=new_param
                    else:
                        out_dict[key]=new_param

        # for multiplication with normal matrices
        else:
            out_dict = {}
            for param_i,const_i in self.param_dict.items():
                new_param = np.matmul(const_i,other)
                out_dict[param_i]=new_param


    def transpose(self):
        transp_self = copy.deepcopy(self)
        for param_i,const_i in transp_self.param_dict.items():
            transp_self.param_dict[param_i]=np.transpose(const_i)
        return transp_self