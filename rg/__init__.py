from . import diagrams, interaction

from . interaction import T, INV_MEASURE, dim, L


#convenient to have these
import operator
import functools
def _product_(l): return functools.reduce(operator.mul, l,1)
def _sum_(terms): return functools.reduce(operator.add, terms , 0)

from . import graph
from .integration import integral
