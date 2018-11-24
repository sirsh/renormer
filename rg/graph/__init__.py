colours = ["blue", "green", "grey", "orange"]
species_labels = ["a", "b", "c", "d"]

from .incidence_matrix import incidence_matrix, _delta_
from .circuit import circuits,set_circuit_flow,set_source_sink_flow,apply_flow
from .viz import show_spanning_trees, ring_diagram, star,arc
from .propagators import scalar_propagator,pole_descriptor
from . import circuit

