import unittest
from cad_kin.roller import Roller
from cad_kin.node import Node
import numpy as np

class Test_Roller(unittest.TestCase):
    def test_constraint(self):
        node1 = Node([0,0])

        nodes = np.array([node1])

        # test 45 degree
        element = {
            "nodes":[0],
            "parametric":False,
            "angle":45
        }
        elem = Roller(element,2)

        constraint = elem(nodes)

        test_vel  = np.array([
            [1],
            [1],
        ])

        self.assertAlmostEquals(0,np.matmul(constraint,test_vel)[0,0])

        test_vel  = np.array([
            [1],
            [0],
        ])

        self.assertNotAlmostEqual(0,np.matmul(constraint,test_vel)[0,0])

        # test 0 degree
        element = {
            "nodes":[0],
            "parametric":False,
            "angle":0
        }
        elem = Roller(element,2)

        constraint = elem(nodes)

        test_vel  = np.array([
            [1],
            [0],
        ])

        self.assertAlmostEqual(0,np.matmul(constraint,test_vel)[0,0])
        
        test_vel  = np.array([
            [1],
            [1],
        ])

        self.assertNotAlmostEqual(0,np.matmul(constraint,test_vel)[0,0])


    def test_constraint_string(self):
        node1 = Node([0,0])

        nodes = np.array([node1])

        element = {
            "nodes":[0],
            "parametric":False,
            "angle":45
        }
        elem = Roller(element,2)

        strings = elem.get_constraint_strings(nodes) 
        # just at edge of contact
        self.assertEqual(
            strings[0],
            '( -70710*v0  +70710*v1 )/100000==0'
        )

    def test_constraint_string_p(self):
        node1 = Node([0,0])

        nodes = np.array([node1])

        element = {
            "nodes":[0],
            "parametric":True,
            "angle":45,
            "p_options":{
                "roll_direction":"any"
            }
        }
        elem = Roller(element,2)

        strings = elem.get_constraint_strings(nodes) 
        # just at edge of contact
        self.assertEqual(
            strings[0],
            '( 70710*v0*a0  +70710*v1*a1 )/100000==0'
        )

        


if __name__ == '__main__':
    unittest.main()