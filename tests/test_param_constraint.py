import unittest
from cad_kin.parametric_constraint import Parameter, ParametricConstraint
import cad_kin.rigid_link as rl
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

    def test2(self):
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


if __name__ == '__main__':
    unittest.main()