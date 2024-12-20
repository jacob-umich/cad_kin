import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import scipy.linalg as sp
from cad_kin.cad_constraint import CadConstraint

class CadTree():
    branch_types = [
        "And",
        "Or"
    ]

    eq_types = [
        "Equal",
        "GreaterEqual",
        "Greater",
        "Less",
        "LessEqual",
        "Inequality"

    ]
    
    def __init__(self,regions,dof,nparam,labels,parameter_info) -> None:
        self.param_info=parameter_info
        self.dof=dof
        self.nparam=nparam
        self.regions = regions
        self.labels = labels
        self.root = CadNode(None,None,0)


    def create_tree(self):
        node_dict = {
            "root":{
                "instances":[self.root],
                "count":1
            }
        }
        # Creates a simple tree that doesn't expand unspecified values or 
        # reorders them
        root = CadNode(None,None,0)
        regions = self.regions
        self.create_nodes(root,regions)
        self.register_nodes(root,node_dict)
        return root

    def register_nodes(self,node,node_dict):
        """Recursively registers all the nodes in a tree in a dictionary"""
        for n in node.children:
            self.register_nodes(n,node_dict)
        if node.id=="root":
            return None
        if node.id in list(node_dict.keys()):
            count = node_dict[node.id]["count"]    
            node.instance = count
            node_dict[node.id]["count"]+=1
            node_dict[node.id]["instances"].append(node)
        else:
            node_dict[node.id]={
                "instances":[node],
                "count":1
            }

    def create_expanded_tree(self):
        regions = self.regions
        labels = self.labels

        node_dict = {
            "root":{
                "instances":[self.root],
                "count":1
            }
        }

        self.create_nodes(self.root,regions)

        # branches trace each path of the tree, creating a separate list of 
        # nodes
        branches = self.scan_branch(self.root,[])

        # check unspecified nodes for each branch. leaf represents the current 
        # progress of building the new tree
        for branch in branches:
            branch_ids = [x.id for x in branch]
            leaves = []
            for n in labels:
                if not (n in branch_ids):
                    # manually create a node
                    if leaves:
                        new_leaves = []
                        while leaves:
                            leaf = leaves.pop()
                            node1 = self.create_man_node(n,"Less")
                            node2 = self.create_man_node(n,"Equal")
                            node3 = self.create_man_node(n,"Greater")
                            leaf.children.append(node1)
                            leaf.children.append(node2)
                            leaf.children.append(node3)
                            node1.parent = leaf
                            node2.parent = leaf
                            node3.parent = leaf
                            new_leaves.append(node1)
                            new_leaves.append(node2)
                            new_leaves.append(node3)
                        leaves=new_leaves

                    else:
                        node1 = self.create_man_node(n,"Less")
                        node2 = self.create_man_node(n,"Equal")
                        node3 = self.create_man_node(n,"Greater")
                        branch[-1].children.append(node1)
                        branch[-1].children.append(node2)
                        branch[-1].children.append(node3)
                        node1.parent = branch[-1]
                        node2.parent = branch[-1]
                        node3.parent = branch[-1]
                        leaves.append(node1)
                        leaves.append(node2)
                        leaves.append(node3)
        # re scan the tree and make new branches to incorporate the newly 
        # created nodes
        branches = self.scan_branch(self.root,[])

        # creating a new tree by re-arranging nodes in each branch based on 
        # specified order of labels, where some nodes exist in multiple 
        # branches. 
        for branch in branches:
            leaf = self.root
            for n in labels:
                for node in branch:
                    if node.id ==n:
                        # check if node is already in new tree by being in a 
                        # previously specified branch. If node has a new 
                        # parent, it has already been recorded in the new tree 
                        # from a previous branch.
                        if node.new_parent==leaf:
                            leaf = node
                            break

                        # check if node is similar to one already on new tree. 
                        # Needs to check heritage too.
                        node = self.check_dup_node(node,leaf,node_dict) # this needs to take in leaf to determine 

                        # check if child is already parent's children
                        if not(node in leaf.new_children):
                            leaf.new_children.append(node)

                        # heritage is the same if already considered node is here, so overwrite is not bad.
                        node.new_parent=leaf
                        leaf = node
                        break
        
        for v in list(self.nodes.values()):
            for node in v["instances"]:
                node.parent = node.new_parent
                node.children = node.new_children

        self.branches = self.scan_branch(self.root,[])

    
    def print_regions(self,fig_path,b_kinematics=False):
        root = self.create_tree()
 
        def delve(node, size_dict):
            width = 0
            depth = 0
            if len(node.children)==0:
                size_dict[f"{node.id}_{node.instance}"]=1
                return size_dict,1,1
            for n in node.children:
                size_dict,width_i, depth_i = delve(n,size_dict)
                width+=width_i
                depth = max(depth_i+1,depth)
            size_dict[f"{node.id}_{node.instance}"]= width
            return size_dict,width,depth
        size_dict = {}
        size_dict,fig_height,fig_width = delve(root,size_dict)
        cell_width = 2.5
        cell_height = 48/72

        fig,ax = plt.subplots(dpi=300,figsize=[(fig_width+2)*cell_width,(fig_height+2)*cell_height]) 
        ax.axis('off')
        ax.margins(0)
        ax.set_xlim(-0.2*cell_width,(fig_width+1)*cell_width)
        ax.set_ylim(-(fig_height+1)*cell_height,cell_height)
        def recurse_print(node,start_pos):
            ax.text(
                start_pos[0],
                cell_height/2+start_pos[1],
                node.__repr__()
            )
            prev_width = 0
            for n in node.children:
                recurse_print(
                    n,
                    [
                        start_pos[0]+cell_width,
                        start_pos[1]-prev_width*cell_height
                    ]
                )
                prev_width+=size_dict[f"{n.id}_{n.instance}"]
        recurse_print(root,[0,0])
        
        fig.savefig(fig_path,format="svg",bbox_inches="tight")

    
    def create_man_node(self,label,eq_type):
        const = CadConstraint(None,self.dof)
        const.eq_type = eq_type
        ind = int(label[1:])
        const.primary=ind
        const.coef = np.array([0])
        const.ind = np.array([-1])
        if "x" in label:
            const.parameter=False
            const.matrix_constraint = np.zeros((1,self.dof))
            const.matrix_constraint[0,ind]=-1
        else:
            const.parameter=True
            
        return CadNode(const,[],0)
    
    def scan_branch(self,node,history):
        histories = []
        history.append(node)
        for child in node.children:
            hist_out = self.scan_branch(child,copy.copy(history))
            if isinstance(hist_out[0],list):
                histories+=hist_out
            else:
                histories.append(hist_out)
        if not node.children:
            return history
        return histories

    def create_nodes(self,parent,region):
        if region.head.name == "Or":
            for r in region:
                self.create_nodes(parent,r)
        if region.head.name =="And":
            current_parent = parent
            for r in region:
                current_parent = self.create_nodes(current_parent,r)
        if region.head.name in self.eq_types:
            new_const = CadConstraint(region,self.dof)
            child = CadNode(new_const,parent,0)
            parent.children.append(child)
            return child
        
    def check_dup_node(self,node,leaf,node_dict):

        if node.id in list(node_dict.keys()):
            count = node_dict[node.id]["count"]
            instances = node_dict[node.id]["instances"]
            matches = [x.const==node.const for x in instances]

            if any(matches):
                # mutliple matches are possible because there is node copying. 
                # need to check heritage of each match
                match_instances = [
                    i for (i, v) 
                    in zip(instances, matches) 
                    if v
                ]
                for match in match_instances:
                    if self.heritage_check(leaf,match):
                        return match
                
                # create new node if node needs to be copied across different 
                # heritages. Needs to copy exact constraint but be a different 
                # node
                if any([node==match for match in match_instances]):
                    node = CadNode(node.const,None,count)
                node.instance = count
                node_dict[node.id]["count"]+=1
                node_dict[node.id]["instances"].append(node)
                return node
            else:
                node.instance = count
                node_dict[node.id]["count"]+=1
                node_dict[node.id]["instances"].append(node)
                return node
        else:
            node_dict[node.id]={
                "instances":[node],
                "count":1
            }
            return node
        
    def heritage_check(self,leaf,match):
        if leaf.id=="root":
            return True
        if leaf.const==match.new_parent.const:
            return self.heritage_check(leaf.new_parent,match.new_parent)
        else:
            return False
        
    def get_constraint_matrix(self,branch):
        matricies = []
        mat = np.concatenate([node.const.matrix_constraint for node in branch[:1:-1] if not(node.const.parameter)],axis=0)
        for node in branch[1:]:
            if ("Less" in node.const.eq_type or "Greater" in node.const.eq_type ) and not(node.const.parameter):
                slack_column = np.zeros((mat.shape[0],1))
                pos = int(node.id[1:])
                if "Less" in node.const.eq_type:
                    slack_column[pos,0]=1
                if "Greater" in node.const.eq_type:
                    slack_column[pos,0]=-1
                mat_i = np.concatenate([mat,slack_column],axis=1)
                matricies.append(mat_i)
            if (node.const.eq_type == "Inequality") and not(node.const.parameter):
                ind = int(node.id[1:])
                new_mats = []
                for mat_i in matricies:
                    mat_i=copy.copy(mat_i)
                    mat_i[ind,:]=np.concatenate([node.const.alt_constraint ,[[0]]],axis=1)
                    new_mats.append(mat_i)
                matricies+=new_mats

        slack_column = np.zeros((mat.shape[0],1))
        mat_i = np.concatenate([mat,slack_column],axis=1)
        matricies.append(mat_i)
        return matricies
    
    def get_velocity(self,mat):
        q,r,p = sp.qr(np.transpose(mat),mode="full",pivoting=True)
        vel = q[:,-1]
        if vel[-1]==0:
            print('here')
        vel = vel/vel[-1]
        vel = vel/np.linalg.norm(vel)
        return vel
    
    def plot_hm(self,ax,label_ax,mat,color,max_v,min_v):
        labels = self.labels[::-1]
        labels = labels[0:-self.nparam]+["b"]
        for i,label in enumerate(labels):
            if "x" in label:
                l_i = label.replace('x',"$v_{")
                l_i+="}$"
                labels[i]=l_i
        im = ax.imshow(mat,vmax=max_v,vmin=min_v)
        ax.set_xticks(np.arange(len(labels)), labels=labels)
        ax.set_yticks([])
        c=plt.Circle((0.5,0.3),0.2,facecolor=color,edgecolor='black',linewidth=2)
        label_ax.add_patch(c)
        label_ax.axes.set_aspect('equal')
        label_ax.axis("off")

        return im
    
    def plot_params(self,signature,ax):
        ax.axis("off")
        
        eqs = signature.get_design()
  
        count=0
        for d in eqs:
            ax.text(0,1-count*0.05-0.05,d)
            count+=1

    def eq_writer(self,nodes,ind,coef,id,symbol): 
        out = "$"
        s = f"{id}"
        s = s[0:1]+"_"+s[1:]+" "
        out+=s
        out+=symbol+" "
        if (ind==-1).all():
            out+= f"{coef[0]}"
        else:
            first=True
            for x,y in zip(coef,ind):
                if y!=-1:
                    s = f"{nodes[y-self.dof].id}"
                    s = s[0:1]+"_"+s[1:]+" "
                    if x>0 and not first:
                        out+="+"
                        first==False
                    if x==1:
                        out+=f"{s}"
                        continue
                    out+= f"{x} {s}"
        out+="$"
        return out
    
    def plot_branch(self,struct,folder_path,branch_mats,branch_n,signature,v_min,v_max):
            vels = []
            for i, mat in enumerate(branch_mats):
                comp=1
                v = self.get_velocity(mat)[0:-1]

                # dont draw redundant velocities
                for v_i in vels:
                    comp = min(np.sum((v-v_i)**2),comp)
                if comp<10**-8:
                    continue
                vels.append(v)
            cell_width = 0.25
            n_mats = len(vels)
            width = self.dof*cell_width*1.05
            n_mat_per_row = int(8*0.975/width)
            cb_width = 0.025*8
            cb_mat_ratio = int(width/cb_width)
            mat_row_gw = cb_mat_ratio
            gw = (cb_mat_ratio*n_mat_per_row+1)*7
            n_struct = 4*n_mat_per_row
            top_mosaic = ["params" for x in range(int(gw/7))]+["struct" for x in range(int(gw/7*6))]
            n_rows = int(np.ceil((n_mats)/n_mat_per_row))
            height = 5+self.dof*cell_width*n_rows*1.2
            n_top = int(20/(self.dof*cell_width*1.2))*5
            mosaic = [top_mosaic]*n_top
            for i in range(n_rows):
                mosaic_i = []
                mosaic_i_label = []
                for j in range(n_mat_per_row):
                    if j+i*n_mat_per_row >= n_mats:
                        mosaic_i_label += [f"." for l in range(cb_mat_ratio*7)]
                        mosaic_i+=[f"." for l in range(cb_mat_ratio*7)]
                    else:
                        mosaic_i_label += [f"hm_label_{j+i*n_mat_per_row}" for l in range(cb_mat_ratio*7)]
                        mosaic_i += [f"hm_{j+i*n_mat_per_row}" for l in range(cb_mat_ratio*7)]
                mosaic_i+=["bar"]*7
                mosaic_i_label+=["bar"]*7
                mosaic+=[mosaic_i_label for l in range(4)]
                mosaic+=[mosaic_i for l in range(16)]

            fig = mpl.figure.Figure(figsize = [8,height],dpi=720)
            ax_dict = fig.subplot_mosaic(mosaic)
            p_nodes = signature.nodes
            p_vals = self.sample_params(p_nodes)
            ims = []
            # max_v = 0
            # min_v = 0
            # for i, mat in enumerate(branch_mats[0:-1]):
            #     max_v = max(np.max(mat),max_v)
            #     min_v = min(np.min(mat),min_v)
            vels = []
            counter = 0
            for mat in branch_mats:
                comp=1  
                v = self.get_velocity(mat)[0:-1]

                # dont draw redundant velocities
                for v_i in vels:
                    comp = min(np.sum((v-v_i)**2),comp)
                if comp<10**-8:
                    continue

                vels.append(v)
                color = np.random.rand(3,)
                struct.move(v*2)
                if sum(v**2)<1e-8:
                    struct.plot(ax_dict["struct"],p_vals,1)
                else:
                    struct.plot(ax_dict["struct"],p_vals,.25,color)
                struct.move(-v*2)
                if sum(v**2)<1e-8:
                    color="white"
                im_i = self.plot_hm(ax_dict[f"hm_{counter}"],ax_dict[f"hm_label_{counter}"],mat,color,v_max,v_min)
                ims.append(im_i)
                counter+=1

            self.plot_params(signature,ax_dict["params"])


            ratio = struct.y_dim/struct.d0
            if ratio<8:
                scale = 8/ratio
                xmin,xmax,ymin,ymax=ax_dict["struct"].axis()
                ax_dict["struct"].axis([xmin*scale,xmax*scale,ymin*scale,ymax*scale])
            fig.colorbar(im_i,ax_dict["bar"])


            fig.savefig(folder_path+f"/branch_{branch_n}.svg", format="svg")
    
    def plot_results(self,struct,folder_path,v_min=-10,v_max=10):
        all_branch_mats = []
        signatures = []

        for k, branch in enumerate(self.branches):
            branch_mats=self.get_constraint_matrix(branch)
            signature = {"included":[]}
            nodes = []
            for node in branch[1:]:
                if node.const.parameter:
                    nodes.append(node)
                else:
                    break
            signature = Signature(nodes,self.dof,self.param_info)
            matched=False
            # going to need to get all branches for one set of params together first
            for i,signature_i in enumerate(signatures):
                if signature==signature_i:
                    all_branch_mats[i]+=branch_mats
                    matched=True

            if not matched:                    
                signatures.append(signature)
                all_branch_mats.append(branch_mats)

        final_branch_mats = []
        final_signatures=[]
        for k,(branch_mats, signature) in enumerate(zip(all_branch_mats,signatures)):
            # check similarity to other results
            matched=False
            for i,signature_i in enumerate(final_signatures):
                # compare matricies
                mask = self.matrix_unique_mask(branch_mats,final_branch_mats[i])
                reverse_mask = self.matrix_unique_mask(final_branch_mats[i],branch_mats)
                # if params are redundantly similar, matricies can be added if they are not similar
                # have to represent all results in this signature
                if signature.deep_equality(signature_i.get_all_sig()):
                    # filter out matrices that are the same
                    branch_mats = [mat for j,mat in enumerate(branch_mats) if mask[j]]

                    final_branch_mats[i]+=branch_mats
                    
                    matched=True
                    signature_i.add_included(signature)

                    break

                # if all matrices are the same while params are different, do not make another fig but keep track of params. mostly here for rigidity cases.
                if not mask.any() and not reverse_mask.any():
                    signature_i.add_included(signature)
                    matched=True
                    break
            
            if not matched:                    
                final_signatures.append(signature)
                final_branch_mats.append(branch_mats)

        path = folder_path+"/results.json"

        out_dict = []
        for m,signature in enumerate(final_signatures):
            out_dict.append(signature.to_json(m))

        with open(path,'w') as f:
            json.dump(out_dict,f, indent=4)

        for k,(branch_mat, signature) in enumerate(zip(final_branch_mats,final_signatures)):
            # try:
                self.plot_branch(struct,folder_path,branch_mat,k,signature,v_min,v_max)
            # except Exception as e:
            #     print(e)

    def get_branch_params(self,branch):
        out = []
        for node in branch[1:]:
            if node.const.parameter:
                out.append(node)
        return out
    
    def matrix_unique_mask(self,mats1,mats2:list):
            comp_mat = np.array(mats2)[:,None,:,:]
            mat = np.array(mats1)
            sums = np.sum((comp_mat-mat)**2,axis=(-1,-2))
            return (sums>=1e-8).any(axis=0)
    
    def signature_equality(self,sig_1,sig_2):
        for i in range(self.nparam):
            node_1 = self.nodes[f"c{i}"]["instances"][sig_1[f"c{i}"]]
            node_2 = self.nodes[f"c{i}"]["instances"][sig_2[f"c{i}"]]
            if node_1.const==node_2.const:
                continue
            if node_1.const.eq_type=="Inequality" or node_2.const.eq_type=="Inequality":
                print("unexpected constraint")
            if self.param_info[i]=="bin":
                if node_1.const.eq_type=="Equal" and node_2.const.eq_type=="Equal":
                    continue
                if (not node_2.const.eq_type=="Equal") and (not node_2.const.eq_type=="Equal"):
                    continue
                # if it gets to here, it is different enough and should have its own plot
                return False
            
            if self.param_info[i]=="flip":
                # if ("Equal" in node_1.const.eq_type) and ("Equal" in node_2.const.eq_type) :
                #     continue

                return False

            if self.param_info[i]=="a":

                return False
            if self.param_info[i]=="b":

                return False
        return True
    
    def sample_params(self,param_nodes):
        params = []
        for p in range(self.nparam):
            node = param_nodes[p]
            if node.const.eq_type=="Inequality":
                if (node.const.ind_l==-1).all():
                    v = node.const.coef_l[0]
                else:
                    v = 0
                    for x,y in zip(node.const.coef_l,node.const.ind_l):
                        if y!=-1:
                            v+=x*params[y-self.dof]
                    v+=1*np.random.random(1)[0]
                # compute other v
                if (node.const.ind_u==-1).all():
                    v2 = node.const.coef_u[0]
                else:
                    v2 = 0
                    for x,y in zip(node.const.coef_u,node.const.ind_u):
                        if y!=-1:
                            v2+=x*params[y-self.dof]
                while v>v2:
                    v+=-.1
                params.append(v)
            
            else:
                if (node.const.ind==-1).all():
                    v = node.const.coef[0]
                else:
                    v = 0
                    for x,y in zip(node.const.coef,node.const.ind):
                        if y!=-1:
                            v+=x*params[y-self.dof]
                if "Less" in node.const.eq_type:
                    v+=-1*np.random.random(1)[0]
                if "Greater" in node.const.eq_type:
                    v+=1*np.random.random(1)[0]
                params.append(v)
        param_vals = np.array(params)
        return param_vals

        

# change to node
class CadNode():

    def __init__(self,const,parent,instance) -> None:
        self.parent = parent
        self.children = []
        self.new_children = []
        self.new_parent = None
        self.instance = instance
        if const:
            self.const = const
            if self.const.parameter:
                out = "a"
            else:
                out = "v"
            out +=f"{self.const.primary}"
            self.id = out
        else:
            self.id = "root"
        
    def __repr__(self) -> str:
        out = self.id
        if getattr(self,"const",False):
            out+= ": "+ f"{self.const}"
        return out
        
    

def get_params(self, branches):
    params = []
    for branch in branches:
        if branch.is_rigid():
            params.append(branch.parameters)
            current_set = np.zeros(self.nparam)
            current_set[0]=1
            for param in range(1,self.nparam):
                found = False
                for param_eq in branch:
                    if param==param_eq.primary:
                        current_set[param]=np.sum(param_coef_i*current_set[param_ind_i] for param_coef_i,param_ind_i in zip(param_eq.coef,param_eq.ind)) 
                        if "Great" in param_eq.eq_type:
                            current_set[param]+=1
                        if "Less" in param_eq.eq_tyep:
                            current_set[param]-=1
                        found = True
                        break
                if not found:
                    current_set[param]="arb"


class Signature():
    symbol_dict = {
        "Less":"<",
        "LessEqual":"\leq",
        "Greater":">",
        "GreaterEqual":"\geq",
        "Equal":"=",
        }
    symbol_dict_basic = {
        "Less":"<",
        "LessEqual":"<=",
        "Greater":">",
        "GreaterEqual":">=",
        "Equal":"=",
        }
    def __init__(self,nodes,dof,param_info):
        self.param_info=param_info
        self.ndof = dof
        self.nodes = nodes
        self.instances = {}
        self.included = []
        for node in self.nodes:
            self.instances[node.id]=node.instance
    
    def get_design(self,basic=False):
        out = []
        if basic:
            symbol_dict=self.symbol_dict_basic
        else:
            symbol_dict=self.symbol_dict

        for node in self.nodes:
            if node.const.eq_type=="Inequality":
                out += [self.gen_equation(
                    self.nodes,
                    node.const.ind_l,
                    node.const.coef_l,
                    node.id,
                    symbol_dict[node.const.l_eq],
                    basic
                )]
            
                out += [self.gen_equation(
                    self.nodes,
                    node.const.ind_u,
                    node.const.coef_u,
                    node.id,
                    symbol_dict[node.const.u_eq],
                    basic
                )]
            
            else:
                out += [self.gen_equation(
                    self.nodes,
                    node.const.ind,
                    node.const.coef,
                    node.id,
                    symbol_dict[node.const.eq_type],
                    basic
                )]
        if basic:
            return "".join(out)
        else:
            return out
    
    def gen_equation(self,nodes,ind,coef,id,symbol,basic="false"):
        if basic:
            out=""
            s = f"{id} "
            out+=s
            out+=symbol+" "
            if (ind==-1).all():
                out+= f"{coef[0]} "
            else:
                first=True
                for x,y in zip(coef,ind):
                    if y!=-1:
                        s = f"{nodes[y-self.ndof].id} "
                        if x>0 and not first:
                            out+="+"
                            first==False
                        if x==1:
                            out+=f"{s}"
                            continue
                        out+= f"{x} {s}"
            return out+"| "
        
        else:
            out = "$"
            s = f"{id}"
            s = s[0:1]+"_"+s[1:]+" "
            out+=s
            out+=symbol+" "
            if (ind==-1).all():
                out+= f"{coef[0]}"
            else:
                first=True
                for x,y in zip(coef,ind):
                    if y!=-1:
                        s = f"{nodes[y-self.ndof].id}"
                        s = s[0:1]+"_"+s[1:]+" "
                        if x>0 and not first:
                            out+="+"
                            first==False
                        if x==1:
                            out+=f"{s}"
                            continue
                        out+= f"{x} {s}"
            out+="$"
            return out

    
    def __eq__(self, other) -> bool:
        return self.instances==other.instances
    
    def deep_equality(self,other):
        if isinstance(other,list):
            truths = []
            for o in other:
                truths.append(self.deep_equality(o))
            return any(truths)
        else:
            for i, (node_1,node_2) in enumerate(zip(self.nodes,other.nodes)):
                if node_1.const==node_2.const:
                    continue
                if node_1.const.eq_type=="Inequality" or node_2.const.eq_type=="Inequality":
                    print("unexpected constraint")
                if self.param_info[i]=="bin":
                    if node_1.const.eq_type=="Equal" and node_2.const.eq_type=="Equal":
                        continue
                    if (not node_1.const.eq_type=="Equal") and (not node_2.const.eq_type=="Equal"):
                        continue
                    # if it gets to here, it is different enough and should have its own plot
                    return False
                
                if self.param_info[i]=="flip":
                    # if ("Equal" in node_1.const.eq_type) and ("Equal" in node_2.const.eq_type) :
                    #     continue

                    return False

                if self.param_info[i]=="a":

                    return False
                if self.param_info[i]=="b":

                    return False
            return True
        
    def get_all_sig(self):
        return [self]+self.included
    
    def add_included(self,other):
        self.included.append(other)
    
    def __repr__(self) -> str:
        return self.get_design(True)
    
    def to_json(self,number=None):
        out = self.instances
        if self.included:
            formulas = []
            for sig in self.included:
                formulas.append(sig.get_design(True))
    
            out["included"] = formulas
        out["values"]=self.get_design(True)
        if number:
            out["id"]=number
        return out