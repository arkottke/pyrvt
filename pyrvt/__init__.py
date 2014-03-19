
from pyrvt import tools
from pyrvt import motions
from pyrvt import peak_calculators

def test(level=1, verbosity=1):
    from numpy.testing import Tester
    return Tester().test(level, verbosity)
