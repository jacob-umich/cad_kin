import matplotlib as mpl
import numpy as np

from cad_kin.structure import Structure
import os

input_path = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
    "structures"
)
output_path = os.path.abspath(os.path.dirname(__file__))

# Four Bar
def four_bar():
    data_path =os.path.join(
            input_path,
            "4_bar.json"
        )
    out_path =os.path.join(
            output_path,
            "example_output",
            "4_bar"
        )
    
    os.makedirs(out_path,exist_ok=True)
    bar = Structure()
    bar.load(data_path)
    tree = bar.cad()
    tree.print_regions(
        os.path.join(
            out_path,
            "design_space.svg",
        )
    )

    # fig = mpl.figure.Figure(dpi=720)
    # ax = fig.add_subplot()
    # bar.plot(ax,np.array([1,1]),1,0.75,annotate=True)
    # fig.savefig(out_path+"/layout.svg", format="svg")

    # tree.plot_results(bar,out_path,-2,2)
    


# Five Bar
def five_bar():
    bar = rigid_kin.Structure('program/rigid_kin/5_bar/5_bar.json')
    bar.gen_constraints()
    bar.add_parametric_contact(1,"any")
    bar.add_parametric_lock(3,4,2)

    fig = mpl.figure.Figure(dpi=720)
    ax = fig.add_subplot()
    bar.plot(ax,np.array([0,0,0,0]),1,annotate=True)
    fig.savefig("program/rigid_kin/5_bar/layout.svg", format="svg")

    labels = bar.get_all_labels()
    geometry = bar.cad('program/rigid_kin/5_bar/5_bar.wls')
    tree = parse_cad.CadTree(geometry,10,4,labels,bar.parameters)
    tree.plot_results(bar,"./program/rigid_kin/5_bar",-3,3)
    bar.session.terminate()

# Six Bar
def six_bar():
    data_path =os.path.join(
            input_path,
            "6_bar.json"
        )
    out_path =os.path.join(
            output_path,
            "example_output",
            "6_bar"
        )
    
    os.makedirs(out_path,exist_ok=True)
    bar = Structure()
    bar.load(data_path)
    # print(bar.compile_constraints())
    tree = bar.cad()
    tree.print_regions(
        os.path.join(
            out_path,
            "design_space.svg",
        )
    )

    # fig = mpl.figure.Figure(dpi=720)
    # ax = fig.add_subplot()
    # bar.plot(ax,np.array([0,0,0,0]),1,hinge_size=0.75)
    # fig.savefig(format="svg")
    # bar.add_parametric_strut(0,3)
    # bar.add_parametric_strut(0,4)
    # bar.add_parametric_lock(4,3,5)
    # bar.add_parametric_contact(1,"any")
    # # bar.add_parametric_contact(2,"any")


# four_bar()
# five_bar() 
six_bar()