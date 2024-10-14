import unittest
from cad_kin.midspan_connect import MidspanConnect
from cad_kin.node import Node
import numpy as np

class Test_Midspan(unittest.TestCase):
    def test_constraint(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,1,2],
            "parametric":False
        }
        elem = MidspanConnect(element,6)

        constraints = elem(nodes)
        constraint1 = constraints[0]
        constraint2 = constraints[1]

        vel  = np.array([
            [0],
            [0],
            [0],
            [1],
            [0],
            [2],
        ])

        self.assert_(np.matmul(constraint2,vel)==0)
        vel  = np.array([
            [0],
            [0],
            [1],
            [0],
            [1],
            [0],
        ])

        self.assert_(np.matmul(constraint1,vel)!=0)

        vel  = np.array([
            [0],
            [0],
            [0],
            [1],
            [0],
            [1],
        ])

        self.assert_(np.matmul(constraint2,vel)!=0)

    def test_constraint_string(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,1,2],
            "parametric":False,
        }
        elem = MidspanConnect(element,6)
        constraints = elem(nodes)
        strings = elem.get_constraint_strings(nodes) 


        # rotation constraint
        cs = [int(abs(constraint)*100000) for constraint in constraints[0]]
        self.assertEqual(
            strings[0],
            f'( -{cs[0]}*v0  +{cs[2]}*v2 )/100000==0'
        )

        # rotation constraint
        cs = [int(abs(constraint)*100000) for constraint in constraints[1]]
        self.assertEqual(
            strings[1],
            f'( {cs[0]}*v0  +{cs[1]}*v1  -{cs[3]}*v3  -{cs[4]}*v4  +{cs[5]}*v5 )/100000==0'
        )


    def test_constraint_string_p(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,1,2],
            "parametric":True,
        }
        elem = MidspanConnect(element,6)
        constraints = elem(nodes)
        strings = elem.get_constraint_strings(nodes) 


        # rotation constraint
        cs = [int(abs(constraint)*200000) for constraint in constraints[0]]
        self.assertEqual(
            strings[0],
            f'( -{cs[0]}*v0*(1-2*a0+a0^2)  +{cs[2]}*v4*(1-2*a0+a0^2) )/100000==0'
        )

        # rotation constraint
        cs = [int(abs(constraint)*100000) for constraint in constraints[1]]
        self.assertEqual(
            strings[1],
            f'( {cs[0]}*v0*(1-3*a0+3*a0^2-a0^3)  +{cs[1]}*v1*(1-3*a0+3*a0^2-a0^3)  -{cs[2]}*v2*(1-3*a0+3*a0^2-a0^3)  -{cs[3]}*v3*(1-3*a0+3*a0^2-a0^3)  -{cs[4]}*v4*(1-3*a0+3*a0^2-a0^3)  +{cs[5]}*v5*(1-3*a0+3*a0^2-a0^3) )/100000==0'
        )

        self.assertEqual(
            strings[2],
            "a0>=0"
        )
        self.assertEqual(
            strings[3],
            "a0<=1"
        )

if __name__ == '__main__':
    unittest.main()