import re
import numpy as np

# This version of the CAD constraint assumes that all design variables are 
# linear with respect to kinematic variables. If for some reason, design variables
# are in the denominator or in higher orders, this program will need to be
# modified.

class CadConstraint():
    operations = [
        "Times",
        "Rational",
        "Plus"
    ]
    symbol_dict = {
        "Equal" :"=",
        "GreaterEqual":">=",
        "Greater":">",
        "Less":"<",
        "LessEqual":"<="
    }
    expr = re.compile("Global`[a-z]([0-9]+)")
    param_expr = re.compile("Global`a([0-9]+)")
    def __init__(self,constraint,total_dof) -> None:
        
        if not constraint:
            self.n_dof = total_dof
            
        else:
            self.n_dof = total_dof
            self.constraint = constraint
            self.eq_type = constraint.head.name
            if self.eq_type=="Inequality":
                self.primary =int(self.expr.search(constraint[2].name).group(1))
                    
                self.coef_l,self.ind_l = self.parse_exp(constraint[0],np.array([1]))
                self.coef_u,self.ind_u = self.parse_exp(constraint[4],np.array([1]))
                self.l_eq = constraint[1].name.replace("Less","Greater")
                self.u_eq = constraint[3].name
                search_res = self.param_expr.search(constraint[2].name)
                if search_res:

                    self.parameter=True
                else:
                    self.parameter=False
                    self.matrix_constraint = np.zeros((1,total_dof))
                    self.alt_constraint = np.zeros((1,total_dof))
                    self.matrix_constraint[0,self.primary]=1
                    self.alt_constraint[0,self.primary]=1
                    if not (self.ind_l==-1).all():
                        self.matrix_constraint[0,self.ind_l]=-self.coef_l
                    if not (self.ind_u==-1).all():
                        self.alt_constraint[0,self.ind_u]=-self.coef_u
            else:
                coef,ind = self.parse_exp(constraint[1],np.array([1]))
                self.coef = coef
                self.ind = ind
            
                # Check if it is a parametric constraint
                search_res = self.param_expr.search(constraint[0].name)
                if search_res:
                    self.primary =int(search_res.group(1))
                    self.parameter=True
                else:
                    self.primary =int(self.expr.search(constraint[0].name).group(1))
                    self.parameter=False

                    self.matrix_constraint = np.zeros((1,total_dof))
                    self.matrix_constraint[0,self.primary]=1
                    if not (ind==-1).all():
                        self.matrix_constraint[0,self.ind]=-coef


    def parse_exp(self,exp,coef):
        if hasattr(exp,"name"):
            if "Global" in exp.name:
                const_ind = int(self.expr.search(exp.name).group(1))
                # add extra indicies if its a parameter
                search_res = self.param_expr.search(exp.name)
                if search_res:
                    return np.array([1]),np.array([const_ind+self.n_dof])
                else:
                    return np.array([1]),np.array([const_ind])
        if hasattr(exp,"head"):

            # -1 for an indice indicates that the coeficient is a scalar
            if exp.head.name == "Rational":
                a,a_ind = self.parse_exp(exp[0],coef)
                b,b_ind = self.parse_exp(exp[1],coef)
                if (a_ind==-1).all() and (b_ind==-1).all():
                    return coef*a/b,a_ind

            if exp.head.name == "Times":
                a,a_ind = self.parse_exp(exp[0],coef)
                b,b_ind = self.parse_exp(exp[1],coef)
                coef=b*a*coef
                if (a_ind==-1).all():
                    return coef,b_ind
                if (b_ind==-1).all():
                    return coef,a_ind
            if exp.head.name == "Plus":
                coefs = []
                inds = []
                for e in exp:
                    c,i = self.parse_exp(e,coef)
                    coefs.append(c)
                    inds.append(i)
                coef = coef*np.concatenate(coefs,axis=0)
                ind = np.concatenate(inds,axis=0)
                return coef,ind
            
            if exp.head.name =="Power":
                print(exp)

        if isinstance(exp,int):
            return coef*np.array([exp]),np.array([-1])
        
    def __eq__(self,other):
        if self.eq_type=="Inequality":
            if not other.eq_type=="Inequality":
                return False

            cond_1 = self.parameter == other.parameter 
            cond_2 = np.linalg.norm(self.coef_l-other.coef_l)<0.001   
            cond_2 = cond_2 and np.linalg.norm(self.coef_u-other.coef_u)<0.001   
            cond_3 = (self.ind_l==other.ind_l).all()
            cond_3 = cond_3 and (self.ind_u==other.ind_u).all()
            cond_4 = self.eq_type==other.eq_type  

        else:
            if other.eq_type=="Inequality":
                return False

            cond_1 = self.parameter == other.parameter 
            cond_2 = np.linalg.norm(self.coef-other.coef)<0.001   
            cond_3 = (self.ind==other.ind).all()
            cond_4 = self.eq_type==other.eq_type  
        return (cond_1 and cond_2 and cond_3 and cond_4) 
    
    def __repr__(self)->str:
        if self.parameter:
            name = f"a{self.primary}"
        else:
            name = f"v{self.primary}"
        if self.eq_type=="Inequality":
            out = self.gen_values("left") + self.symbol_dict["LessEqual"]
            out += name
            out += self.gen_values("right") + self.symbol_dict["LessEqual"]

        else:
            out = name +self.symbol_dict[self.eq_type]
            out += self.gen_values()
        return out
    
    def gen_values(self,coef_type="normal"):
        match coef_type:
            case "normal":
                coef = self.coef
                ind = self.ind
            case "left":
                coef = self.coef_l
                ind = self.ind_l
            case "right":
                coef = self.coef_u
                ind = self.ind_u
        
        
        out = ""
        if (ind==-1).all():
            out+= f"{coef[0]}"

        else:
            first=True
            for x,y in zip(coef,ind):
                if y!=-1:
                    if y>=self.n_dof:
                        var_type = "a"
                        shift = self.n_dof
                    else:
                        var_type = "v"
                        shift = 0
                    s = f"{var_type}{y-shift}"
                    if x>0 and not first:
                        out+="+"
                    if x==1:
                        out+=f"{s}"
                        first=False
                        continue
                    first=False
                    out+= f"{x:.2f} {s}"
        return out
