import unittest
import cad_kin.parametric_constraint
import numpy as np

class Test_Strut(unittest.TestCase):
    def test1(self):
        p1 = cad_kin.parametric_constraint.Parameter("a1")
        p2 = cad_kin.parametric_constraint.Parameter("a2")
        p3 = cad_kin.parametric_constraint.Parameter("a3")
        c = cad_kin.parametric_constraint.Parameter("1")
        self.assert_((c*p2).symbol=="a2")
        self.assert_((p2*c).symbol=="a2")

        self.assert_((p1*p2).symbol=="a1*a2")
        self.assert_((p1*p2*p3).symbol=="a1*a2*a3")
        print((p1*p2*p3*p1).symbol)
        self.assert_((p1*p2*p3*p1).symbol=="a1*a1*a2*a3")


if __name__ == '__main__':
    unittest.main()