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

import logging as log
import numpy.random as npr
import madsenlab.ctmixtures.analysis as analysis
from base_rule import BaseInteractionRule
from collections import defaultdict
import random
import numpy.random as npr

class ConformistCopyingRule(BaseInteractionRule):
    """
    Implements a neutral copying process via Moran dynamics, taking an instance of a lattice model at construction.
    Returns control to the caller after each step(), so that other code can run to determine completion,
    take samples, etc.
    """

    def __init__(self, model):
        self.model = model
        self.sc = self.model.simconfig

    def step(self, agent, timestep):
        """
        Implements a single time step in the neutral drift Moran model, starting from a focal agent,
        and then one of the focal agent's neighbors at random (this rule knows nothing about
        how "neighbors" are represented, so the rule itself is fully generic to many
        population structures, including those with long-distance connections.

        """

        prob = 0

        if npr.random() < self.sc.conformism_strength:
            # execute a local conformism rule among neighbors
                # choose a random locus
            num_loci = self.sc.num_features
            rand_locus = random.randint(0,num_loci)

            # get the traits from all neighbors at that locus
            neighbor_ids = self.model.get_all_neighbors_for_agent(agent)
            trait_cnts = defaultdict(int)
            for neighbor in neighbor_ids:
                trait_cnts[self.model.agentgraph[neighbor]["agent"].traits[rand_locus]] += 1
            sorted_cnts = sorted(trait_cnts, reverse=True)
            log.debug("sorted traits: %s", sorted_cnts)

            # the most frequent trait will be the first item in the sorted trait list
            selected_trait = sorted_cnts[0]
            log.debug("selected trait: %s", selected_trait)
            agent.traits[rand_locus] = selected_trait

        else:
            # execute a normal random copy
            neighbor = self.model.get_random_neighbor_for_agent(agent.id)
            differing_features = analysis.get_different_feature_positions_locusallele(agent.traits, neighbor.traits)
            old_agent_traits = list(agent.traits)
            if len(differing_features) == 1:
                random_feature = differing_features[0]
            else:
                rand_feature_num = npr.randint(0, len(differing_features))
                random_feature = differing_features[rand_feature_num]
            neighbor_trait = neighbor.traits[random_feature]
            agent.traits[random_feature] = neighbor_trait
            log.debug("agent %s: old: %s  neighbor: %s  post: %s differing: %s feature: %s val: %s ", agent.id, old_agent_traits, neighbor.traits, agent.traits,differing_features, random_feature, neighbor_trait )

        # track the interaction and time
        self.model.update_interactions(timestep)


class AntiConformistCopyingRule(BaseInteractionRule):
    """
    Implements a neutral copying process via Moran dynamics, taking an instance of a lattice model at construction.
    Returns control to the caller after each step(), so that other code can run to determine completion,
    take samples, etc.
    """

    def __init__(self, model):
        self.model = model
        self.sc = self.model.simconfig

    def step(self, agent, timestep):
        """
        Implements a single time step in the neutral drift Moran model, starting from a focal agent,
        and then one of the focal agent's neighbors at random (this rule knows nothing about
        how "neighbors" are represented, so the rule itself is fully generic to many
        population structures, including those with long-distance connections.

        """

        if npr.random() < self.sc.anticonformism_strength:
            # execute a local conformism rule among neighbors
                # choose a random locus
            num_loci = self.sc.num_features
            rand_locus = random.randint(0,num_loci)

            # get the traits from all neighbors at that locus
            neighbor_ids = self.model.get_all_neighbors_for_agent(agent)
            trait_cnts = defaultdict(int)
            for neighbor in neighbor_ids:
                trait_cnts[self.model.agentgraph[neighbor]["agent"].traits[rand_locus]] += 1
            sorted_cnts = sorted(trait_cnts, reverse=False)
            log.debug("sorted traits: %s", sorted_cnts)

            # the least frequent trait will be the first item in the sorted trait list
            selected_trait = sorted_cnts[0]
            log.debug("selected trait: %s", selected_trait)
            agent.traits[rand_locus] = selected_trait

        else:
            # TODO:  Reevaluate "different positions" code for both conformist and neutral - remove it.  Only appropriate for axelrod models

            # execute a normal random copy
            neighbor = self.model.get_random_neighbor_for_agent(agent.id)
            differing_features = analysis.get_different_feature_positions_locusallele(agent.traits, neighbor.traits)
            old_agent_traits = list(agent.traits)
            if len(differing_features) == 1:
                random_feature = differing_features[0]
            else:
                rand_feature_num = npr.randint(0, len(differing_features))
                random_feature = differing_features[rand_feature_num]
            neighbor_trait = neighbor.traits[random_feature]
            agent.traits[random_feature] = neighbor_trait
            log.debug("agent %s: old: %s  neighbor: %s  post: %s differing: %s feature: %s val: %s ", agent.id, old_agent_traits, neighbor.traits, agent.traits,differing_features, random_feature, neighbor_trait )

        # track the interaction and time
        self.model.update_interactions(timestep)





