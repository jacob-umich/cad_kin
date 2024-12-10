from cad_kin.roller import Roller
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

class ContactBC(Roller):
    eq_symbol = ">="
    def __init__(self,element,n_dof) -> None:
        super().__init__(element,n_dof)

    def __call__(self,nodes):
        return super().__call__(nodes)

    def plot(self,nodes,drawing_thickness,drawing_color ='#D0D0D0',params=None ):
        if self.b_parametric:
            if (params==None).all():
                return []
            if (params==0).all():
                return []
            self.angle = super().get_angle_from_params(params)
            if params[1]<0:
                self.angle+=180
        pos, dofs = self.get_node_info(nodes)
        x = pos[0]
        y = pos[1]
        t = drawing_thickness*0.99
        roller = []

        base = plt.Rectangle((x-0.7*t,y-t*.5),width=t*1.4,height=t*0.1,facecolor='#D0D0D0')
        roller.append(base)

        for i in range(5):
            start = (i-2)*0.3-0.2
            mark = plt.Rectangle((x+start*t,y-t*0.8),width=t*.42,height=t*0.06,angle = 45, facecolor='#D0D0D0')
            roller.append(mark)

        for s in roller:
            transf = mpl.transforms.Affine2D().rotate_deg_around(x,y,self.angle)
            s.set_transform(transf)
            if self.b_parametric:
                s.set_facecolor("#389ac7ff")
        return roller
    