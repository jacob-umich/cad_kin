import unittest
from cad_kin.parametric_constraint import Parameter, ParametricConstraint
import cad_kin.rigid_link as rl
from cad_kin import Pin, Structure
from cad_kin.node import Node
from cad_kin.linear_parametric_node import LinearParametricNode

import numpy as np

class Test_Strut(unittest.TestCase):
    def test1(self):
        p1 = Parameter("a1")
        p2 = Parameter("a2")
        p3 = Parameter("a3")
        c = Parameter("1")
        self.assert_((c*p2).symbol=="a2")
        self.assert_((p2*c).symbol=="a2")

        self.assert_((p1*p2).symbol=="a1*a2")
        self.assert_((p1*p2*p3).symbol=="a1*a2*a3")
        print((p1*p2*p3*p1).symbol)
        self.assert_((p1*p2*p3*p1).symbol=="a1*a1*a2*a3")

    def testlinkCall_parametricNode(self):
        node1 = Node([0,0])
        node3 = Node([1,-1])
        node4 = Node([1,1])
        node2 = LinearParametricNode(node3,node4,"a0")
    

        nodes = np.array([node1,node2,node3,node4])

        element = {
            "nodes":[0,1],
            "parametric":False,
        }
        elem = rl.RigidLink(element,6)

        constraint = elem(nodes)

    def testpinCall(self):
        node1 = Node([0,0])
    

        nodes = np.array([node1])

        element = {
            "nodes":[0],
            "parametric":True,
        }
        elem = Pin(element,2)

        constraint = elem(nodes)

    def testsimple_struct(self):
        in_dict = {
            "nodes": [
                [
                    0.0,
                    0.0
                ],
                [
                    0.0,
                    5.0
                ],
                [
                    5.0,
                    0.0
                ],
                [
                    5.0,
                    5.0
                ],
                [
                    2.5,
                    0.0
                ]
            ],
            "elements": [
                {
                    "type":"link",
                    "nodes":[0,2],
                    "thickness":0.75
                },
                {
                    "type":"link",
                    "nodes":[1,3],
                    "thickness":0.75
                },       
                {
                    "type":"link",
                    "nodes":[3,4],
                    "thickness":0.75
                },       
                {
                    "type":"pin",
                    "nodes":[0]
                },
                {
                    "type":"roller",
                    "angle":90,
                    "nodes":[1]
                },
                {
                    "type":"midspan",
                    "nodes":[0,4,2],
                    "parametric":True
                }
            ]
        }

    
        st = Structure(in_dict)

        print(st.compile_constraints())


if __name__ == '__main__':
    unittest.main()