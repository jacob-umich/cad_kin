import rigid_kin
import parse_cad
import matplotlib as mpl
import numpy as np
def six_bar():

    bar = rigid_kin.Structure('program/rigid_kin/6_bar/6_bar.json')
    bar.gen_constraints()
    # bar.add_parametric_strut(0,3)
    # bar.add_parametric_strut(0,4)
    # bar.add_parametric_lock(3,2,4)
    # bar.add_parametric_lock(4,3,5)
    # bar.add_parametric_contact(1,"any")
    # bar.add_parametric_contact(2,"any")

    labels = bar.get_all_labels()
    geometry = bar.cad('program/rigid_kin/6_bar/6_bar.wls')
    bar.session.terminate()

def two_bar():
    bar = rigid_kin.Structure('program/rigid_kin/2_bar/2_bar.json')
    bar.gen_constraints()
    bar.add_parametric_roller(2,"y")
    bar.add_parametric_strut(0,2)
    bar.add_parametric_contact(1,"any")
    bar.add_parametric_contact(2,"any")
    bar.add_parametric_contact(2,"any")
    bar.add_parametric_contact(2,"any")
    bar.add_parametric_lock(1,2,0)
    bar.add_parametric_roller(1,"y")

    labels = bar.get_all_labels()
    geometry = bar.cad('program/rigid_kin/2_bar/2_bar.wls')
    bar.session.terminate()

def rigid():
    bar = rigid_kin.Structure('program/rigid_kin/rigid_triangles/t_8.json')
    bar.gen_constraints()
    geometry = bar.cad('program/rigid_kin/rigid_triangles/out.wls')
    bar.session.terminate()

def results():
    six_bar_times = [1.625041,1.820198,2.106280,3.026143,5.120045,6.1,14.8,328.493634,69108.6]
    two_bar_times = [1.542401,1.542436,1.558767,1.591277,1.639045,2.017968,2.899364,6.647081,15.630903] #3875.994638 for 18
    rigid_times = [1.673596,1.692530,1.729237,1.586826,1.665556,1.673303,1.760703,1.742394]
    n_design_params = [0,1,2,3,4,5,6,7,8]
    n_rigid = [2*x+6 for x in range(8)]
    six_params = [x+12 for x in n_design_params]
    two_params = [x+6 for x in n_design_params]
    fig = mpl.figure.Figure(dpi=720)
    ax = fig.add_subplot()
    ax.plot(n_design_params,six_bar_times,label="six-bar",color="#005a00ff",linewidth=4)
    ax.plot(n_design_params,two_bar_times,label="two-bar",color="#389ac7ff",linewidth=4)
    ax.set_yscale('log')
    ax.legend()
    fig.savefig("program/rigid_kin/scaling.svg", format="svg")

    fig = mpl.figure.Figure(dpi=720)
    ax = fig.add_subplot()
    ax.plot(six_params,six_bar_times,label="six-bar",color="#005a00ff",linewidth=4)
    ax.plot(two_params,two_bar_times,label="two-bar",color="#389ac7ff",linewidth=4)
    ax.plot(n_rigid,rigid_times,label="rigid",color="#e5472eff",linewidth=4)
    ax.set_yscale('log')
    ax.legend()
    fig.savefig("program/rigid_kin/scaling_params.svg", format="svg")

two_bar()