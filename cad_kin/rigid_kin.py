import numpy as np
import json
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wlexpr
from dotenv import load_dotenv
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PatchCollection
import matplotlib.transforms as transforms
import os
import logging

load_dotenv()



def nullprojector(a,b):
    a_unit = a/np.linalg.norm(a)
    a_unit = a_unit[:,None]
    b = b[:,None]
    proj = b-np.matmul(np.matmul(a_unit,np.transpose(a_unit)),b)
    return proj

def check_inline(a,b):
    a_unit = a/np.linalg.norm(a)
    a_unit = a_unit[:,None]
    b = b[:,None]
    p = np.matmul(np.matmul(a_unit,np.transpose(a_unit)),b)
    return np.linalg.norm(p)<np.linalg.norm(a) and np.matmul(np.transpose(p),a_unit)>0

def mid_link_condition(node_pos,elem_pos,d0):
    a = elem_pos[2:4]-elem_pos[:2]
    b = node_pos-elem_pos[:2]
    c = node_pos-elem_pos[2:4]

    if np.linalg.norm(b)==0 or  np.linalg.norm(c)==0:
        return False
    b_orth_a = nullprojector(a,b)
    dist = np.linalg.norm(b_orth_a)
    cond = abs(dist-d0)<1e-3 and check_inline(a,b)
    if cond:
        return True
    else:
        return False
    
def mid_link_contact(node_pos,elem_pos,d0,n,node_dof,elem_dof):
    a = elem_pos[2:4]-elem_pos[:2]
    b = node_pos-elem_pos[:2]
    da = np.zeros((2,n))
    db = np.zeros((2,n))
    for i in range (2):
        db[i,node_dof[i]]=1
        db[i,elem_dof[i]]=-1
        da[i,elem_dof[i+2]]=1
        da[i,elem_dof[i]]=-1
    p = nullprojector(a,b) 
    q = a/np.linalg.norm(a)
    p = q*np.dot(q,b)  
    dq = (np.linalg.norm(a)*da-np.matmul(a[:,None],np.matmul(a[None,:],da))/np.linalg.norm(a))/np.linalg.norm(a)**2
    dp = dq*np.dot(q,b)+np.matmul(q[:,None],(np.matmul(b[None,:],dq)+np.matmul(q[None,:],db)))
    constraint = np.matmul((b-p)[None,:],(db-dp))
    return constraint

def node_to_node_condition(node_pos,elem_pos,d0):

    d = node_pos-elem_pos
    cond = abs(np.linalg.norm(d)-d0)<1e-3 and np.linalg.norm(d)>0
    return cond
    
def node_to_node_contact(node_pos,elem_pos,n,node_dof,elem_dof):
    a = node_pos-elem_pos
    da = np.zeros((2,n))
    for i in range (2):
        da[i,node_dof[i]]=1
        da[i,elem_dof[i]]=-1
    return np.matmul(a[None,:],da)/np.linalg.norm(a)

def rot_lock(pos,dofs,ndof):
    a = pos[2:4]-pos[0:2]
    b = pos[4:6]-pos[0:2]
    da = np.zeros((2,ndof))
    db = np.zeros((2,ndof))
    for i in range (2):
        db[i,dofs[i+4]]=1
        db[i,dofs[i]]=-1
        da[i,dofs[i+2]]=1
        da[i,dofs[i]]=-1
    a_norm = np.linalg.norm(a)
    b_norm = np.linalg.norm(b)
    a_ip = a_norm**2
    a_t = a[None,:]
    a = a[:,None]
    b_t = b[None,:]
    b = b[:,None]

    const = np.matmul(a_t,db)/(a_norm*b_norm)-np.matmul(a_t,b)*np.matmul(b_t,db)/(a_norm*b_norm**3)+np.matmul(b_t,da)/(b_norm*a_norm)-np.matmul(b_t,a)*np.matmul(a_t,da)/(b_norm*a_norm**3)
    return const


def mid_link_connect(pos,dofs,ndof):
    ref = np.array([[3/5],[4/5]])
    a = pos[2:4]-pos[0:2]
    b = pos[4:6]-pos[0:2]
    da = np.zeros((2,ndof))
    db = np.zeros((2,ndof))
    for i in range (2):
        db[i,dofs[i+4]]=1
        db[i,dofs[i]]=-1
        da[i,dofs[i+2]]=1
        da[i,dofs[i]]=-1
    a_norm = np.linalg.norm(a)
    b_norm = np.linalg.norm(b)
    a_ip = a_norm**2
    a_t = a[None,:]
    a = a[:,None]
    b_t = b[None,:]
    b = b[:,None]
    const_1 = ((np.matmul(b_t,da)+np.matmul(a_t,db))*a_norm-np.matmul(a_t,b)*np.matmul(a_t,da)/(a_norm))/(a_ip)
    const_2 = np.matmul(np.transpose(ref),
        (
            db/b_norm
        )-(
            np.matmul(b,np.matmul(b_t,db))/b_norm**3
        )-(
            da/a_norm
        )+(
            np.matmul(b,np.matmul(b_t,db))/b_norm**3
        )
    )
    if np.all(const_2==0):
        ref = np.array([[4/5],[3/5]])
        const_2 = np.matmul(np.transpose(ref),
            (
                db/b_norm
            )-(
                np.matmul(b,np.matmul(b_t,db))/b_norm**3
            )-(
                da/a_norm
            )+(
                np.matmul(b,np.matmul(b_t,db))/b_norm**3
            )
        )
    return const_1, const_2


class Structure():
    def __init__(self, input_file):
        with open(input_file,'r') as f:
            data = json.load(f)
        self.nodes_posns = np.array(data['nodes'])
        nodes = []
        coord_number = 0
        for node_i in data['nodes']:
            node_i_out = {
                'pos':np.array(node_i,dtype='float64'),
                'dofs':np.array([coord_number,coord_number+1])
            }
            coord_number+=2
            nodes.append(node_i_out)
        self.nodes = np.array(nodes)
        self.nodes_dofs = np.array([x['dofs'] for x in self.nodes])
        self.n_dof = self.nodes_dofs.size
        links = []
        for el in data['elems']:
            pos = np.concatenate([self.nodes_posns[el[0]],self.nodes_posns[el[1]]],axis=0)
            dofs = np.concatenate([self.nodes_dofs[el[0]],self.nodes_dofs[el[1]]],axis=0)
            
            link_i = {
                'pos':pos,
                'dofs':dofs
            }
            links.append(link_i)
        self.links = links

        self.straights = data.get('straights',{})
        self.bc=data['boundary_conds']
        self.d0=data['thickness']
        self.p_constraint = []
        self.parameters = []
        self.n_params = 0

        # get plotting info
        maxes = np.max(self.nodes_posns,axis=0)
        self.x_dim =maxes[0]
        self.y_dim =maxes[1]
        self.apsect = self.x_dim/self.y_dim
        self.scale = 0.1/self.d0

        try:
            self.session = WolframLanguageSession(os.getenv("WOLFRAM_KERNEL_PATH"),kernel_loglevel=logging.DEBUG)
        except Exception:
            print('wolfram kernel not initialized')
    
    def gen_constraints(self):
        link_consts = []
        for i in self.links:
            pos = i['pos']
            constr_i = np.zeros((1,self.n_dof))
            constr_i[0,i['dofs']]=np.array([
                [
                    pos[0]-pos[2],
                    pos[1]-pos[3],
                    -(pos[0]-pos[2]),
                    -(pos[1]-pos[3]),
                ],
            ])
            link_consts.append(constr_i)
        straight_const = []
        for i in self.straights:
            node_pos = self.nodes[i.get("node")].get("pos")
            elem_pos = self.links[i.get("element")].get("pos")
            pos = np.concatenate([elem_pos,node_pos],axis=0)
            node_dofs = self.nodes[i.get("node")].get("dofs")
            elem_dofs = self.links[i.get("element")].get("dofs")
            dof = np.concatenate([elem_dofs,node_dofs],axis=0)
            const_1, const_2 = mid_link_connect(pos,dof,self.n_dof)
            straight_const.append(const_1)
            straight_const.append(const_2)
        bc_const = []
        for i in self.bc:
            const_i = np.zeros((1,self.n_dof))
            const_i[0,i]=1
            bc_const.append(const_i)
    
        contact_const = []
        for node in self.nodes:
            for elem in self.links:
                const = self.gen_contact_constraint(node,elem)
                if const.any():
                    contact_const.append(const)

        link_consts=np.concatenate(link_consts,axis=0)
        bc_const = np.concatenate(bc_const,axis=0)
        if contact_const:
            contact_const = np.concatenate(contact_const,axis=0)

        self.contact_constraints=contact_const
        if straight_const:
            straight_const = np.concatenate(straight_const,axis=0)
            self.constraints = np.concatenate([link_consts,straight_const,bc_const], axis=0)
        else:
            self.constraints = np.concatenate([link_consts,bc_const], axis=0)

    def gen_contact_constraint(self,node,elem):

                if mid_link_condition(node['pos'],elem['pos'],self.d0):

                    return mid_link_contact(
                            node['pos'],
                            elem['pos'],
                            self.d0,
                            self.n_dof,
                            node['dofs'],
                            elem['dofs']
                        )
                else:
                    for i in range(2):
                        if node_to_node_condition(node['pos'],elem['pos'][i*2:i*2+2],self.d0):
                            return node_to_node_contact(node['pos'],elem['pos'][i*2:i*2+2],self.n_dof,node['dofs'],elem['dofs'][i*2:i*2+2])
                        else:
                            out =  np.zeros((1,self.n_dof))
                    return out
    
    def cad(self,output_file,precision=16):
        out = "out=CylindricalDecomposition[\n{"

        for constraint in self.constraints:
            first=True
            out+="("
            for i,v in enumerate(constraint):
                if v!=0:
                    if v<0:
                        s = "-"
                    else:
                        s = "+"
                    if( first and s=="+"):
                        out+=f"  {int(abs(v)*10**(precision))}*x{i}"
                        first=False
                    else:
                        out+=f"  {s}{int(abs(v)*10**(precision))}*x{i}"
                        first = False
            out+=f")/{10**precision}==0,\n"

        for constraint in self.contact_constraints:
            out+="("
            first=True
            for i,v in enumerate(constraint):
                if v!=0:
                    if v<0:
                        s = "-"
                    else:
                        s = "+"
                    if( first and s=="+"):
                        out+=f"  {int(abs(v)*10**(precision))}*x{i}"
                        first=False
                    else:
                        out+=f"  {s}{int(abs(v)*10**(precision))}*x{i}"
                        first = False
            out+=f")/{10**precision}>=0,\n"


        parameters = []
        for k,p_c_i in enumerate(self.p_constraint):
            out+="("
            map_i = p_c_i["map"]
            const_i = p_c_i["const"]
            first=True
            for p in map_i.keys():
                if p in parameters:
                    continue
                else: 
                    parameters.append(p)

            for i,v in enumerate(const_i[0]):
                if v!=0:
                    param_string = ""
                    for p,v_ind in map_i.items():
                        if i in v_ind:
                            param_string+=p
                    if v<0:
                        s = "-"
                    else:
                        s = "+"
                    if( first and s=="+"):
                        out+=f"  {int(abs(v)*10**(precision))}*x{i}*{param_string}"
                        first=False
                    else:
                        out+=f"  {s}{int(abs(v)*10**(precision))}*x{i}*{param_string}"
                        first = False
            if self.parameters[k]["type"]=="strut" or self.parameters[k]["type"]=="contact":
                out+=f")/{10**precision}>=0,\n"   
            else:    
                out+=f")/{10**precision}==0,\n"   


        out=out[0:-2]+out[-1:]
        out +="},\n{"

        for k in parameters:
            out+=f"{k}, "
        for i in range(self.n_dof):
            if i!=self.n_dof-1:
                out+=f"x{self.n_dof-1-i}, "
            else:
                out+=f"x{self.n_dof-1-i}"+"}\n"
        out+=']'
        with open(output_file,"w") as f:
            f.write(out)
        try:
            return self.session.evaluate(wlexpr(out))
        except Exception as e:
            print(e)
            print('wolfram kernel not initialized')

    def get_dof_labels(self):
        out = []
        for i in range(self.n_dof):
            out.append(f"x{self.n_dof-1-i}")
        return out
    
    def get_all_labels(self,):
        dofs = self.get_dof_labels()
        params = []
        for i in range(self.n_params):
            params.append(f"c{i}")
        return params+dofs
    
    def plot_hinge(self,point):
        x = point[0]
        y = point[1]
        t = self.d0*0.99
        hinge = []
        r = plt.Rectangle((x-0.7*t,y-t*0.6),width=t*1.4,height=t*0.6,facecolor='#D0D0D0')
        hinge.append(r)
        c = plt.Circle((x,y),t*.7,facecolor='#D0D0D0') 
        base = plt.Rectangle((x-0.9*t,y-t*0.7),width=t*1.8,height=t*0.1,facecolor='#D0D0D0')
        hinge.append(c)
        hinge.append(base)

        for i in range(5):
            start = (i-2)*0.3-0.2
            mark = plt.Rectangle((x+start*t,y-t),width=t*.42,height=t*0.06,angle = 45, facecolor='#D0D0D0')
            hinge.append(mark)
        return hinge
    
    def plot_roller(self,point):
        x = point[0]
        y = point[1]
        t = self.d0*0.99
        roller = []
        r = plt.Rectangle((x-0.7*t,y-t*0.7),width=t*1.4,height=t*0.7,facecolor='#D0D0D0')
        roller.append(r)
        c = plt.Circle((x,y),t*.7,facecolor='#D0D0D0')
        
        roll_1 = plt.Circle((x-0.3*t,y-t*0.9),t*.2,facecolor='#D0D0D0')

        roll_2 = plt.Circle((x+t*0.3,y-t*0.9),t*.2,facecolor='#D0D0D0') 

        base = plt.Rectangle((x-0.7*t,y-t*1.2),width=t*1.4,height=t*0.1,facecolor='#D0D0D0')
        roller.append(c)
        roller.append(base)
        roller.append(roll_1)
        roller.append(roll_2)

        for i in range(5):
            start = (i-2)*0.3-0.2
            mark = plt.Rectangle((x+start*t,y-t*1.5),width=t*.42,height=t*0.06,angle = 45, facecolor='#D0D0D0')
            roller.append(mark)
        return roller

    def plot_contact(self,point):
        x = point[0]
        y = point[1]
        t = self.d0*0.99
        roller = []

        base = plt.Rectangle((x-0.7*t,y-t*.5),width=t*1.4,height=t*0.1,facecolor='#D0D0D0')
        roller.append(base)

        for i in range(5):
            start = (i-2)*0.3-0.2
            mark = plt.Rectangle((x+start*t,y-t*0.8),width=t*.42,height=t*0.06,angle = 45, facecolor='#D0D0D0')
            roller.append(mark)
        return roller

    def plot_link(self,elem,t,color):
        x2=elem['pos'][0]
        x1=elem['pos'][2]
        y2=elem['pos'][1]
        y1=elem['pos'][3]
        if x1==x2:
            angle=90*np.sign(y2-y1)
        else:
            slope = (y2-y1)/(x2-x1)
            angle = np.arctan(slope)*180/np.pi
            if (x2-x1)<0:
                angle+=180

        length = np.sqrt((y2-y1)**2+(x2-x1)**2)
        r = plt.Rectangle((x1,y1-t/2),width=length,height=t,angle=angle,rotation_point=(x1,y1),facecolor=color,edgecolor=color,linewidth=2)
        return r
    
    def plot_strut(self,elem,t,color,cable=False):
        x2=elem['pos'][0]
        x1=elem['pos'][2]
        y2=elem['pos'][1]
        y1=elem['pos'][3]
        if x1==x2:
            angle=90*np.sign(y2-y1)
        else:
            slope = (y2-y1)/(x2-x1)
            angle = np.arctan(slope)*180/np.pi
            if (x2-x1)<0:
                angle+=180

        length = np.sqrt((y2-y1)**2+(x2-x1)**2)
        if cable:
            # r=mpl.lines.Line2D([x1,x2],[y1,y2], linewidth=2,linestyle="--")
            r = plt.Rectangle((x1,y1-t/2),width=length,height=t/10,angle=angle,rotation_point=(x1,y1),facecolor="none",edgecolor=color,linewidth=2)
        else:
            r = plt.Rectangle((x1,y1-t/2),width=length,height=t,angle=angle,rotation_point=(x1,y1),facecolor="none",edgecolor=color,linewidth=2)
        return r
    
    def plot_params(self,param_vals,alpha):
        patches = []


        for param_i in self.parameters:
            if param_i["type"] =="roller":

                val_i = param_vals[param_i["params"]]
                if (val_i==0).all():
                    continue
                else:
                    x = self.nodes_posns[param_i["node"]][0]
                    y = self.nodes_posns[param_i["node"]][1]
                    if len(val_i)==2:
                        if val_i[1]==0:
                            angle=90
                        else:
                            angle = np.arctan(-val_i[0]/val_i[1])*180/np.pi
                    else:
                        if param_i["dir"]=="x":
                            angle = 0
                        else:
                            angle =90

                    roller = self.plot_roller(self.nodes_posns[param_i["node"]])
                    for s in roller:
                        transf = mpl.transforms.Affine2D().rotate_deg_around(x,y,angle)
                        s.set_transform(transf)
                        # ax.add_patch(s)
                    patches+=roller
                    p_ind = param_i["params"][0]
                    # ax.text(x,y+self.d0,f"$c_{p_ind}$",color='cyan')
                    continue
            if param_i["type"] =="contact":

                val_i = param_vals[param_i["params"]]
                if (val_i==0).all():
                    continue
                else:
                    x = self.nodes_posns[param_i["node"]][0]
                    y = self.nodes_posns[param_i["node"]][1]
                    if len(val_i)==2:
                        if val_i[1]==0:
                            angle=90
                        else:
                            angle = np.arctan(-val_i[0]/val_i[1])*180/np.pi
                    else:
                        if param_i["dir"]=="x":
                            angle = 0
                        else:
                            angle =90
                    if val_i[1]<0:
                        angle+=180
                    roller = self.plot_contact(self.nodes_posns[param_i["node"]])
                    for s in roller:
                        transf = mpl.transforms.Affine2D().rotate_deg_around(x,y,angle)
                        s.set_transform(transf)
                        # ax.add_patch(s)
                    patches+=roller
                    p_ind = param_i["params"][0]
                    # ax.text(x,y+self.d0,f"$c_{p_ind}$",color='cyan')
                    continue
            if param_i["type"]=="pin":
                val_i = param_vals[param_i["params"]]
                if val_i==0:
                    continue
                else:
                    x = self.nodes_posns[param_i["node"]][0]
                    y = self.nodes_posns[param_i["node"]][1]
                    pin = self.plot_hinge([x,y])
                    patches+=pin
                    continue

            if param_i["type"]=="link":
                val_i = param_vals[param_i["params"]]
                if val_i==0:
                    continue
                else:
                    elem = {
                        "pos":np.concatenate([self.nodes[param_i["node_i"]]["pos"],self.nodes[param_i["node_j"]]["pos"]],axis=0)
                    }
                    patches.append(self.plot_link(elem,self.d0*.99,"green"))
                    continue
            if param_i["type"]=="strut":
                val_i = param_vals[param_i["params"]]
                if val_i==0:
                    continue
                else:
                    elem = {
                        "pos":np.concatenate([self.nodes[param_i["node_i"]]["pos"],self.nodes[param_i["node_j"]]["pos"]],axis=0)
                    }
                    if val_i<0:
                        cable=True
                    else:
                        cable=False
                    patches.append(self.plot_strut(elem,self.d0*.99,"green",cable))
                    continue

        param_collection = PatchCollection(patches,match_original=True)
        param_collection.set_facecolor("#389ac7ff")
        param_collection.set_alpha(alpha)
           
        return param_collection
        # param_rollers = filter(lambda x:x["type"]=="roller",self.parameters)
        # for roller in param_rollers:
        #     node = roller["node"]

    def plot_bc(self,alpha):

        patches = []
        t = self.d0*0.99
        for i in range(0,self.n_dof,2):
            if (i in self.bc) and (i+1 in self.bc):
                hinge = self.plot_hinge(self.nodes[int(i/2)]["pos"])
                patches+=hinge
                continue


            if (i in self.bc):
                roller = self.plot_roller(self.nodes[int(i/2)]["pos"])
                for s in roller:
                    transf = mpl.transforms.Affine2D().rotate_deg_around(self.nodes[int(i/2)]["pos"][0],self.nodes[int(i/2)]["pos"][1],90)
                    s.set_transform(transf)
                patches+=roller
                
                continue
            if (i+1 in self.bc):
                roller = self.plot_roller(self.nodes[int(i/2)]["pos"])
                patches+=roller
                continue
        param_collection = PatchCollection(patches,match_original=True)
        param_collection.set_alpha(alpha)
        return param_collection
    
    def plot_struct(self, alpha):
        patches = []
        t = self.d0*0.85
        for elem in self.links:
            patches.append(self.plot_link(elem,t,"#005a00ff"))

        bars = PatchCollection(patches,match_original=True)
        patches = []

        for n in self.nodes_posns:
            c=plt.Circle((n[0],n[1]),t/2,facecolor='white',edgecolor='black',alpha=alpha,linewidth=2)
            patches.append(c)

        nodes = PatchCollection(patches,match_original=True)
        return nodes,bars
    
    def plot(self,ax,param_vals,alpha,color=None,annotate=False):

        # scale = min(6.4/bar.x_dim,4/bar.y_dim)
        # print(scale)            
        
        ax.axis('off')
        ax.set_xlim(self.x_dim*-0.2,self.x_dim*1.3)
        ax.set_ylim(self.y_dim*-0.2,self.y_dim*1.3)
        ax.axes.set_aspect('equal')

        c1 = self.plot_bc(alpha)
        if isinstance(param_vals,np.ndarray):

            c2 = self.plot_params(param_vals,alpha)
        nodes,bars = self.plot_struct(alpha)
        if isinstance(color,np.ndarray):
            c1.set_facecolor(color)
            c1.set_alpha(alpha)
            if isinstance(param_vals,np.ndarray):
                c2.set_facecolor(color)
                c2.set_alpha(alpha)
            bars.set_color(color)       
            bars.set_alpha(alpha)
            nodes.set_facecolor(color)       
            # nodes.set_alpha(alpha)
        
        ax.add_collection(bars)
        if isinstance(param_vals,np.ndarray):
            c2.set_edgecolor(None)
            ax.add_collection(c2)
        
        c1.set_edgecolor(None)

        ax.add_collection(c1)
        ax.add_collection(nodes)
        if annotate:
            self.plot_annotation(ax,param_vals)
        # t = ax.transData.inverted()+transforms.Affine2D().scale(scale)+transforms.Affine2D().translate(1.6,4)+fig.dpi_scale_trans
        # c1.set_transform(c1.get_transform()+t)
    
    def plot_annotation(self,ax,param_vals):
        t = self.d0*.99

        for i, (d,p) in enumerate(zip(self.nodes_dofs,self.nodes_posns)):
            ax.text(p[0],p[1]+t,f"($v_{{{d[0]}}}$, $v_{{{d[1]}}}$)")
            ax.text(p[0]+t,p[1],f"$node_{i}$")
        for i, pval in enumerate(param_vals):
            ax.text(1,1,f"$\\alpha_{i}$",color="cyan")
        


    #     return a_collection
    
    def move(self,positions):
        for i,n in enumerate(self.nodes):
            n['pos']+=positions[n['dofs']]
        self.nodes_posns = np.array([x['pos'] for x in self.nodes])

        for link in self.links:
            link['pos']+=positions[link['dofs']]

    def add_constraint(self,const):
        self.constraints=np.concatenate([self.constraints,const],axis=0)

    def add_parametric_constraint(self,p_map,p_const):
        self.p_constraint.append(
            {
                "map":p_map,
                "const":p_const
            }
        )

    def add_parametric_roller(self,node,roll_direction):
        if roll_direction=="any":
            param_map = {
                f"c{self.n_params}":[self.nodes[node]["dofs"][0]],
                f"c{self.n_params+1}":[self.nodes[node]["dofs"][1]]
            }
            param_const = np.zeros((1,self.n_dof))
            param_const[0,self.nodes[node]["dofs"][0]]=1
            param_const[0,self.nodes[node]["dofs"][1]]=1
            params = [self.n_params,self.n_params+1]
            self.n_params+=2
            param_rule = ["a","b"]

        if roll_direction=="x":
            param_map = {
                f"c{self.n_params}":[self.nodes[node]["dofs"][1]],
            }
            param_const = np.zeros((1,self.n_dof))
            param_const[0,self.nodes[node]["dofs"][1]]=1

            params = [self.n_params]
            self.n_params+=1
            param_rule = ["bin"]

        if roll_direction=="y":
            param_map = {
                f"c{self.n_params}":[self.nodes[node]["dofs"][0]],
            }
            param_const = np.zeros((1,self.n_dof))
            param_const[0,self.nodes[node]["dofs"][0]]=1

            params = [self.n_params]
            self.n_params+=1
            param_rule = ["bin"]

        self.parameters.append({
            "type":"roller",
            "node":node,
            "dir":roll_direction,
            "params":params,
            "rule":param_rule
        })
        
        self.add_parametric_constraint(param_map,param_const)

    def add_parametric_pin(self,node):

        param_map = {
            f"c{self.n_params}":[self.nodes[node]["dofs"][1]],
        }
        param_const = np.zeros((1,self.n_dof))
        param_const[0,self.nodes[node]["dofs"][1]]=1

        self.add_parametric_constraint(param_map,param_const)

        param_map = {
            f"c{self.n_params}":[self.nodes[node]["dofs"][0]],
        }
        param_const = np.zeros((1,self.n_dof))
        param_const[0,self.nodes[node]["dofs"][0]]=1
        params = [self.n_params]
        self.n_params+=1
        
        self.add_parametric_constraint(param_map,param_const)
        
        # additional rules that help with post processing
        self.parameters.append({
            "type":"pin",
            "node":node,
            "params":params,
            "rule":["bin"]
        })

    def add_parametric_link(self,i,j):
        dofs =np.concatenate([self.nodes_dofs[i],self.nodes_dofs[j]],axis=0)
        diff = self.nodes_posns[j]-self.nodes_posns[i]
        p_const = np.zeros((1,self.n_dof))
        p_const[0,[dofs]] = [-diff[0],-diff[1],diff[0],diff[1]]
        param_map = {f"c{self.n_params}":dofs}
        params = [self.n_params]
        self.n_params+=1
        self.add_parametric_constraint(param_map,p_const)
        self.parameters.append({
            "type":"link",
            "node_i":i,
            "node_j":j,
            "params":params,
            "rule":["bin"]
        })
        
    def add_parametric_lock(self,hinge,na,nb):
        self.parameters.append({
            "type":"lock",
            "node_i":hinge,
            "node_j":na,
            "node_k":nb,
            "rule":["bin"]
        })
        # dofs =np.concatenate([self.nodes_dofs[i],self.nodes_dofs[j]],axis=0)
        p_const = self.get_rotation_constraint(hinge,na,nb)

        d1 = self.nodes[hinge].get("dofs")
        d2 = self.nodes[na].get("dofs")
        d3 = self.nodes[nb].get("dofs")
        dofs = np.concatenate([d1,d2,d3],axis=0)

        param_map = {f"c{self.n_params}":dofs}
        self.n_params+=1
        self.add_parametric_constraint(param_map,p_const)

        
    def add_parametric_strut(self,i,j):
        dofs =np.concatenate([self.nodes_dofs[i],self.nodes_dofs[j]],axis=0)
        diff = self.nodes_posns[j]-self.nodes_posns[i]
        p_const = np.zeros((1,self.n_dof))
        p_const[0,[dofs]] = [-diff[0],-diff[1],diff[0],diff[1]]
        param_map = {f"c{self.n_params}":dofs}
        params = [self.n_params]
        self.n_params+=1
        self.add_parametric_constraint(param_map,p_const)
        self.parameters.append({
            "type":"strut",
            "node_i":i,
            "node_j":j,
            "params":params,
            "rule":["flip"]
        })

    def add_parametric_contact(self,node,roll_direction):
        self.add_parametric_roller(node,roll_direction)
        self.parameters[-1]["type"]="contact"

    def get_rotation_constraint(self,nc,na,nb):
        p1 = self.nodes[nc].get("pos")
        p2 = self.nodes[na].get("pos")
        p3 = self.nodes[nb].get("pos")
        d1 = self.nodes[nc].get("dofs")
        d2 = self.nodes[na].get("dofs")
        d3 = self.nodes[nb].get("dofs")
        pos = np.concatenate([p1,p2,p3],axis=0)
        dofs = np.concatenate([d1,d2,d3],axis=0)

        c = rot_lock(pos,dofs,self.n_dof)
        return c
        

    def simple_rigid(self):
        rank = np.linalg.matrix_rank(self.constraints)
        if rank==self.n_dof:
            print("structure is rigid")
            return 0
        q,r = np.linalg.qr(np.transpose(self.constraints),"complete")
        return(q[:,rank:])


