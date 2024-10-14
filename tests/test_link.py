import unittest
from cad_kin.rigid_link import RigidLink
from cad_kin.node import Node
import numpy as np

class Test_Link(unittest.TestCase):
    def test_contact_detection(self):
        node1 = Node([0,0])
        node2 = Node([1,.1])
        node3 = Node([2,0])

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,2],
            "parametric":False,
            "thickness":0.1
        }
        elem = RigidLink(element,6)

        # just at edge of contact
        self.assert_(elem.detect_contact(nodes)[0])

        # not in contact
        node2.pos = np.array([1,1])
        self.assertFalse(elem.detect_contact(nodes)[0])

        # in line with linkage but not in contact
        node2.pos = np.array([2.5,0])
        self.assertFalse(elem.detect_contact(nodes)[0])

        # interpenetration contact
        node2.pos = np.array([1,0.05])
        self.assert_(elem.detect_contact(nodes)[0])

        # should raise error because admissible direction cant be determined explicitly
        node2.pos = np.array([1,0.05])
        try:
            elem.detect_contact(nodes)[0]
            self.assert_(False)
        except:
            self.assert_(True)

    def test_constraint(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,1],
            "parametric":False,
        }
        elem = RigidLink(element,6)

        constraint = elem(nodes)

        vel  = np.array([
            [0],
            [0],
            [0],
            [1],
            [0],
            [2],
        ])

        self.assert_(np.matmul(constraint,vel)==0)

        vel  = np.array([
            [0],
            [0],
            [1],
            [0],
            [1],
            [0],
        ])

        self.assert_(np.matmul(constraint,vel)!=0)

        vel  = np.array([
            [0],
            [0],
            [0],
            [1],
            [0],
            [1],
        ])

        self.assert_(np.matmul(constraint,vel)==0)

    def test_constraint_string(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,1],
            "parametric":False,
        }
        elem = RigidLink(element,6)

        strings = elem.get_constraint_strings(nodes) 
        # just at edge of contact
        self.assertEqual(
            strings[0],
            '( -100000*v0  +100000*v2 )/100000==0'
        )
    def test_constraint_string_p(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,1],
            "parametric":True,
        }
        elem = RigidLink(element,6)

        strings = elem.get_constraint_strings(nodes) 
        # just at edge of contact
        self.assertEqual(
            strings[0],
            '( -100000*v0*a0  +100000*v2*a0 )/100000==0'
        )


if __name__ == '__main__':
    unittest.main()