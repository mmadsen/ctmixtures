#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

import networkx as nx
from numpy.random import RandomState
import logging as log

class LocusAlleleTraitFactory(object):
    """
    A trait factory for models where agents have F loci and T possible traits per locus.
    Individuals are initialized with a list of F random integers, each chosen from 0 to T-1.
    The result is given as a Python list, and stored as the individual's initial trait set.

    This factory may be dynamically loaded from its fully qualified name in a configuration file,
     and passed the simulation configuration object in its constructor.  The instantiating
     code then calls initialize_population(graph), passing it a NetworkX graph of nodes, previously
     constructed
    """

    def __init__(self, simconfig):
        self.simconfig = simconfig
        self.prng = RandomState()  # allow the library to choose a seed via OS specific mechanism

    def initialize_population(self,graph):
        nf = self.simconfig.num_features
        nt = self.simconfig.num_traits
        for nodename in graph.nodes():
            graph.node[nodename]['traits'] = self.prng.randint(0, nt, size=nf)


