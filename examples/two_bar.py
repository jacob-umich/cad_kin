import matplotlib as mpl
import numpy as np

from cad_kin.structure import Structure
import os

input_path = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
    "structures"
)
output_path = os.path.abspath(os.path.dirname(__file__))

data_path =os.path.join(
        input_path,
        "2_bar.json"
    )
out_path =os.path.join(
        output_path,
        "example_output",
        "2_bar"
    )

os.makedirs(out_path,exist_ok=True)
bar = Structure()
bar.load(data_path)
print(bar.compile_constraints(True))
print(bar.get_modes())

fig = mpl.figure.Figure(dpi=720)
ax = fig.add_subplot()
bar.plot(ax,np.array([1,1]),1,0.75,annotate=True)
fig.savefig(out_path+"/layout.svg", format="svg")