import unittest

def grad(d,e):
    return (potential_energy(d+e)-potential_energy(d))/e

class TestGrad(unittest.TestCase):

    def test_grad(self):
        e = Ang*pow(10,-10)
        ## The values of e that verify the test are Ang*pow(10,-8) and Ang*pow(10,-9)
        ## If e is too small we approach a divide by zero error case and hence values no longer match 
        for d in x:
            self.assertAlmostEqual(force(d),-grad(d,e), places = 7)

unittest.main(argv=[''], verbosity=2, exit=False)
