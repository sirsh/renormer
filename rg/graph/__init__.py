colours = ["blue", "green"]
species_labels = ["a", "b"]

from .incidence_matrix import incidence_matrix
from .circuit import circuits,set_circuit_flow,set_source_sink_flow,apply_flow
from .viz import show_spanning_trees, ring_diagram
from .propagators import scalar_propagator,pole_descriptor
from . import circuit

