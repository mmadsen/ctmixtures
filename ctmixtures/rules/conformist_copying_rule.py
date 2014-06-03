#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
This rule implements the original Axelrod model on a lattice, given descriptions in Axelrod (1997) and:

@book{Barrat2009,
    Author = {Barrat, A and Barth\'elemy, M and Vespignani, A},
    Publisher = {Cambridge University Press},
    Title = {Dynamical processes on complex networks},
    Year = {2009}}


"""

from collections import defaultdict
import logging as log
from collections import OrderedDict
import numpy.random as npr
from operator import itemgetter
from base_rule import BaseInteractionRule


class BaseNeighborConformismRule(BaseInteractionRule):
    """
    Implements the common logic for conformist and anticonformist transmission.  If the CONFORMISM_FLAG is
    set to True, then this is pro-conformism, if False, then anticonformism.

    """

    def step(self, agent, timestep):
        """
        Implements a single time step in the neutral drift Moran model, starting from a focal agent,
        and then one of the focal agent's neighbors at random (this rule knows nothing about
        how "neighbors" are represented, so the rule itself is fully generic to many
        population structures, including those with long-distance connections.

        """

        # self.strength is either conformism or anticonformism strength, depending upon which class it is.
        # this allows asymmetric strengths

        if npr.random() < self.strength:
            # execute a local conformism rule among neighbors
                # choose a random locus
            num_loci = self.sc.num_features
            rand_locus = npr.randint(0,num_loci)
            #log.debug("(anti)conformism - random locus: %s with type: %s", rand_locus, self.ruletype)

            # get the traits from all neighbors at that locus
            neighbors = self.model.get_all_neighbors_for_agent(agent.id)
            trait_cnts = defaultdict(int)
            for neighbor in neighbors:
                ntrait = neighbor.traits[rand_locus]
                trait_cnts[ntrait] += 1

            ordered_cnts = sorted(trait_cnts.items(), key=itemgetter(1), reverse=self.CONFORMISM_FLAG)
            target_trait = ordered_cnts[0]
            target_trait_key = target_trait[0]
            target_trait_count = target_trait[1]

            #log.debug("sorted traits: %s", ordered_cnts)
            #log.debug("selected trait: %s", target_trait_key)

            agent.traits[rand_locus] = target_trait_key

        else:
            # execute a normal random copy
            neighbor = self.model.get_random_neighbor_for_agent(agent.id)
            num_loci = self.sc.num_features
            rand_locus = npr.randint(0,num_loci)
            #log.debug("a/conformism but below rate, copy randomly - random locus: %s", rand_locus)

            neighbor_trait = neighbor.traits[rand_locus]
            agent.traits[rand_locus] = neighbor_trait

        # track the interaction and time
        self.model.update_interactions(rand_locus, timestep)






class ConformistCopyingRule(BaseNeighborConformismRule):
    """
    Implements a neutral copying process via Moran dynamics, taking an instance of a lattice model at construction.
    Returns control to the caller after each step(), so that other code can run to determine completion,
    take samples, etc.
    """

    def __init__(self, model):
        self.model = model
        self.sc = self.model.simconfig
        self.CONFORMISM_FLAG = True
        self.ruletype = 'conformism'
        self.strength = self.sc.conformism_strength
        log.debug("Conformist rule operating at strength: %s", self.strength)





class AntiConformistCopyingRule(BaseNeighborConformismRule):
    """
    Implements a neutral copying process via Moran dynamics, taking an instance of a lattice model at construction.
    Returns control to the caller after each step(), so that other code can run to determine completion,
    take samples, etc.
    """

    def __init__(self, model):
        self.model = model
        self.sc = self.model.simconfig
        self.CONFORMISM_FLAG = False
        self.ruletype = 'anticonformism'
        self.strength = self.sc.anticonformism_strength
        log.debug("Anticonformist rule operating at strength: %s", self.strength)




