import unittest
from cad_kin.rotation_lock import RotationLock
from cad_kin.node import Node
import numpy as np

class Test_RotationLock(unittest.TestCase):
    
    # basic rotation tests
    def test_constraint(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,1,2],
            "parametric":False,
        }
        elem = RotationLock(element,6)

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
            [0],
            [1],
            [0],
            [1],
        ])

        self.assert_(np.matmul(constraint,vel)!=0)

    # advanced rotation tests
    def test_constraint2(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([np.sqrt(2)/2,np.sqrt(2)/2])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,1,2],
            "parametric":False,
        }
        elem = RotationLock(element,6)

        constraint = elem(nodes)

        vel  = np.array([
            [0],
            [0],
            [0],
            [1],
            [-np.sqrt(2)/2],
            [np.sqrt(2)/2],
        ])

        self.assert_(np.matmul(constraint,vel)==0)

        vel  = np.array([
            [0],
            [0],
            [0],
            [0],
            [np.sqrt(2)/2],
            [-np.sqrt(2)/2],
        ])

        self.assert_(np.matmul(constraint,vel)!=0)

        # this constraint doesnt control linkage extension
        vel  = np.array([
            [0],
            [0],
            [0],
            [0],
            [np.sqrt(2)/2],
            [np.sqrt(2)/2],
        ])

        self.assertAlmostEqual(np.matmul(constraint,vel)[0,0],0)

    def test_constraint_string(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0,1,2],
            "parametric":False,
        }
        elem = RotationLock(element,6)

        strings = elem.get_constraint_strings(nodes) 

        self.assertEqual(
            strings[0],
            '( 424264*v0  +141421*v1  -282842*v2  -282842*v3  -141421*v4  +141421*v5 )/100000==0'
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
        elem = RotationLock(element,6)

        strings = elem.get_constraint_strings(nodes) 
        # just at edge of contact
        self.assertEqual(
            strings[0],
            '( 424264*v0*a0  +141421*v1*a0  -282842*v2*a0  -282842*v3*a0  -141421*v4*a0  +141421*v5*a0 )/100000==0'
        )


if __name__ == '__main__':
    unittest.main()