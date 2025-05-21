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

    def __repr__(self):
        return self.symbol

    def __hash__(self):
        return hash(self.symbol)
    def __eq__(self, value):
        if isinstance(value,Parameter):
            return self.ids == value.ids
        else:
            return self.symbol==value

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
    precision=10
    def __init__(self, param_dict):
        self.param_dict = param_dict

    def __matmul__ (self, other):
        b_pcosntraint = isinstance(other,ParametricConstraint)

        if b_pcosntraint:
            out_dict = {}
            for param_i,const_i in self.param_dict.items():
                for param_j,const_j in other.param_dict.items():
                    key = param_i*param_j
                    values = out_dict.get(key,[])
                    new_param = np.matmul(const_i,const_j)
                    if len(values)>0:
                        out_dict[key]+=new_param
                    else:
                        out_dict[key]=new_param

        # for multiplication with normal matrices
        else:
            out_dict = {}
            for param_i,const_i in self.param_dict.items():
                new_param = np.matmul(const_i,other)
                out_dict[param_i]=new_param
        return ParametricConstraint(out_dict)
    
    def __mul__(self,other):
        b_pcosntraint = isinstance(other,ParametricConstraint)

        if b_pcosntraint:
            out_dict = {}
            for param_i,const_i in self.param_dict.items():
                for param_j,const_j in other.param_dict.items():
                    key = param_i*param_j
                    values = out_dict.get(key,[])
                    new_param = const_i*const_j
                    if len(values)>0:
                        out_dict[key]+=new_param
                    else:
                        out_dict[key]=new_param

        # for multiplication with normal numbers
        else:
            out_dict = {}
            for param_i,const_i in self.param_dict.items():
                new_param = const_i*other
                out_dict[param_i]=new_param
        return ParametricConstraint(out_dict)      

    def __add__ (self, other):
        b_pcosntraint = isinstance(other,ParametricConstraint)

        if b_pcosntraint:
            out_dict = {}

            # add mutual params
            for param_i,const_i in self.param_dict.items():
                for param_j,const_j in other.param_dict.items():
                    if param_i==param_j:
                        out_dict[param_i]=const_j+const_i

            # get params that are not mutual

            key_i = list(self.param_dict.keys())
            key_j = list(self.param_dict.keys())

            for ki in key_i:
                if ki in key_j:
                    continue
                else:
                    out_dict[ki]=self.param_dict[ki]
            for kj in key_j:
                if kj in key_i:
                    continue
                else:
                    out_dict[kj]=other.param_dict[kj]
            return ParametricConstraint(out_dict)
        # for addition with normal matrices
        else:
            values = self.param_dict.get("1",[])
            if len(values)>0:
                values+=other
            else:
                values=other

            self.param_dict["1"]=values
            
            return self
    def __sub__ (self, other):
        b_pcosntraint = isinstance(other,ParametricConstraint)

        if b_pcosntraint:
            out_dict = {}

            # add mutual params
            for param_i,const_i in self.param_dict.items():
                for param_j,const_j in other.param_dict.items():
                    if param_i==param_j:
                        out_dict[param_i]=const_i-const_j

            # get params that are not mutual

            key_i = list(self.param_dict.keys())
            key_j = list(self.param_dict.keys())

            for ki in key_i:
                if ki in key_j:
                    continue
                else:
                    out_dict[ki]=self.param_dict[ki]
            for kj in key_j:
                if kj in key_i:
                    continue
                else:
                    out_dict[kj]=-other.param_dict[kj]
            return ParametricConstraint(out_dict)
        # for addition with normal matrices
        else:
            values = self.param_dict.get("1",[])
            if len(values)>0:
                values-=other
            else:
                values=-other

            self.param_dict["1"]=values
            
            return self

    def __neg__ (self):

        out_dict = {}

        for param_i,const_i in self.param_dict.items():
            out_dict[param_i]=-const_i

        return ParametricConstraint(out_dict)

    def __getitem__(self, key):
        out_dict = {}
        for k,v in self.param_dict.items():
            out_dict[k]=v[key]
        return ParametricConstraint(out_dict)

    def __setitem__(self, key, value):
        for k,v in self.param_dict.items():
            v[key]=value

    def __delitem__(self, key):
        for k,v in self.param_dict.items():
            del v[key]

    def transpose(self):
        transp_self = copy.deepcopy(self)
        for param_i,const_i in transp_self.param_dict.items():
            transp_self.param_dict[param_i]=np.transpose(const_i)
        return transp_self
    
    def dot(self,other):
        b_pcosntraint = isinstance(other,ParametricConstraint)
        if b_pcosntraint:
            return self@other.transpose()
        else:
            return self@np.transpose(other)
    
    def concatenate(self,other,axis):
        b_pcosntraint = isinstance(other,ParametricConstraint)

        if b_pcosntraint:
            out_dict = {}

            # add mutual params
            for param_i,const_i in self.param_dict.items():
                for param_j,const_j in other.param_dict.items():
                    if param_i==param_j:
                        out_dict[param_i]=np.concatenate([const_i,const_j],axis=axis)

            # get params that are not mutual

            key_i = list(self.param_dict.keys())
            key_j = list(self.param_dict.keys())

            for ki in key_i:
                if ki in key_j:
                    continue
                else:
                    zero_vec = np.zeros(other.param_dict[key_j[0]].shape)
                    out_dict[ki]=np.concatenate([self.param_dict[ki],zero_vec],axis=axis)
            for kj in key_j:
                if kj in key_i:
                    continue
                else:                    
                    zero_vec = np.zeros(other.param_dict[key_i[0]].shape)
                    out_dict[kj]=np.concatenate([zero_vec,self.param_dict[kj]],axis=axis)
            return ParametricConstraint(out_dict)
        # for addition with normal matrices
        else:
            out_dict ={}
            values = self.param_dict.get("1",[])
            values = np.concatenate([values,other],axis=axis)
            out_dict["1"]=values
            zero_vec = np.zeros(shape=other.shape)
            for param_i,const_i in self.param_dict.items():
                if param_i!="1":
                    out_dict[param_i]=np.concatenate([const_i,zero_vec],axis=axis)
            return ParametricConstraint(out_dict)

    def sqrt(self):
        out_dict = {}
        for k,v in self.param_dict.items():
            out_dict[k]=np.sqrt(v)
        return ParametricConstraint(out_dict)
    
    def get_string(self):
        n_constraints = self.param_dict[list(self.param_dict.keys())[0]].shape[0]
        n_dof = self.param_dict[list(self.param_dict.keys())[0]].shape[1]
        strings = []
        for j in range(n_constraints):
            terms = [""]*n_dof
            for k,v in self.param_dict.items():
                if not any(v[j]):
                    continue
                for i,c in enumerate(v[j]):
                    if abs(c)>1e-15:
                        param_string = str(k)
                        if c<0:
                            s = "-"
                        else:
                            s = ""
                        if len(terms[i])>0:
                            terms[i]+=f"+({s}{int(abs(c)*10**(self.precision))})*v{i}*{param_string} "
                        else:
                            terms[i]+=f"({s}{int(abs(c)*10**(self.precision))})*v{i}*{param_string} "

            out = "+".join([t for t in terms if t])
            strings.append(out)
        return strings
    def __str__(self):
        strings = self.get_string()
        return ",\n".join(strings)
    def __repr__(self):
        return str(self)