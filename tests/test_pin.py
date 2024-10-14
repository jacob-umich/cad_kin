import unittest
from cad_kin.pin import Pin
from cad_kin.node import Node
import numpy as np

class Test_Pin(unittest.TestCase):
    
    def test_constraint(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0],
            "parametric":False,
        }
        elem = Pin(element,6)

        constraint = elem(nodes)

        vel  = np.array([
            [0],
            [0],
            [0],
            [1],
            [0],
            [2],
        ])

        self.assert_(np.matmul(constraint,vel)[0]==0)
        self.assert_(np.matmul(constraint,vel)[1]==0)

        vel  = np.array([
            [0],
            [1],
            [1],
            [0],
            [1],
            [0],
        ])

        self.assert_(np.matmul(constraint,vel)[0]==0)
        self.assert_(np.matmul(constraint,vel)[1]!=0)

        vel  = np.array([
            [1],
            [0],
            [0],
            [1],
            [0],
            [1],
        ])

        self.assert_(np.matmul(constraint,vel)[0]!=0)
        self.assert_(np.matmul(constraint,vel)[1]==0)

    def test_constraint_string(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0],
            "parametric":False,
        }
        elem = Pin(element,6)

        strings = elem.get_constraint_strings(nodes) 
        
        self.assertEqual(
            strings[0],
            '( 100000*v0 )/100000==0'
        )
        self.assertEqual(
            strings[1],
            '( 100000*v1 )/100000==0'
        )

    def test_constraint_string_p(self):
        node1 = Node([0,0])
        node2 = Node([1,0])
        node3 = Node([2,0])
    

        nodes = np.array([node1,node2,node3])

        element = {
            "nodes":[0],
            "parametric":True,
        }
        elem = Pin(element,6)

        strings = elem.get_constraint_strings(nodes) 

        self.assertEqual(
            strings[0],
            '( 100000*v0*a0 )/100000==0'
        )
        self.assertEqual(
            strings[1],
            '( 100000*v1*a0 )/100000==0'
        )


if __name__ == '__main__':
    unittest.main()