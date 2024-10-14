import unittest
from cad_kin.contact_boundary import ContactBC
from cad_kin.node import Node
import numpy as np

class Test_ContactBC(unittest.TestCase):
    def test_constraint(self):
        node1 = Node([0,0])

        nodes = np.array([node1])

        # test 45 degree
        element = {
            "nodes":[0],
            "parametric":False,
            "angle":45
        }
        elem = ContactBC(element,2)

        constraint = elem(nodes)

        test_vel  = np.array([
            [1],
            [1],
        ])

        self.assertAlmostEquals(0,np.matmul(constraint,test_vel)[0,0])

        test_vel  = np.array([
            [0],
            [1],
        ])

        self.assertGreaterEqual(np.matmul(constraint,test_vel)[0,0],0)

        test_vel  = np.array([
            [1],
            [0],
        ])

        self.assertLessEqual(np.matmul(constraint,test_vel)[0,0],0)

        # test 0 degree
        element = {
            "nodes":[0],
            "parametric":False,
            "angle":0
        }
        elem = ContactBC(element,2)

        constraint = elem(nodes)

        test_vel  = np.array([
            [1],
            [1],
        ])

        self.assertGreaterEqual(np.matmul(constraint,test_vel)[0,0],0)

        test_vel  = np.array([
            [0],
            [1],
        ])

        self.assertGreaterEqual(np.matmul(constraint,test_vel)[0,0],0)

        test_vel  = np.array([
            [1],
            [0],
        ])

        self.assertAlmostEqual(np.matmul(constraint,test_vel)[0,0],0)

        test_vel  = np.array([
            [0],
            [-1],
        ])

        self.assertLessEqual(np.matmul(constraint,test_vel)[0,0],0)


    def test_constraint_string(self):
        node1 = Node([0,0])

        nodes = np.array([node1])

        element = {
            "nodes":[0],
            "parametric":False,
            "angle":45
        }
        elem = ContactBC(element,2)

        strings = elem.get_constraint_strings(nodes) 
        # just at edge of contact
        self.assertEqual(
            strings[0],
            '( -70710*v0  +70710*v1 )/100000>=0'
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
        elem = ContactBC(element,2)

        strings = elem.get_constraint_strings(nodes) 
        # just at edge of contact
        self.assertEqual(
            strings[0],
            '( 70710*v0*a0  +70710*v1*a1 )/100000>=0'
        )

        element = {
            "nodes":[0],
            "parametric":True,
            "angle":45,
            "p_options":{
                "roll_direction":"x"
            }
        }
        elem = ContactBC(element,2)

        strings = elem.get_constraint_strings(nodes) 
        # just at edge of contact
        self.assertEqual(
            strings[0],
            '( 100000*v1*a0 )/100000>=0'
        )

        element = {
            "nodes":[0],
            "parametric":True,
            "angle":45,
            "p_options":{
                "roll_direction":"y"
            }
        }
        elem = ContactBC(element,2)

        strings = elem.get_constraint_strings(nodes) 
        # just at edge of contact
        self.assertEqual(
            strings[0],
            '( 100000*v0*a0 )/100000>=0'
        )

        


if __name__ == '__main__':
    unittest.main()