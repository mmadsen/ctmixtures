#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
.. module:: simulation_data
    :platform: Unix, Windows
    :synopsis: Data object for storing metadata and parameter about a simulation run in MongoDB, via the Ming ORM.

.. moduleauthor:: Mark E. Madsen <mark@madsenlab.org>

"""

import logging as log
from ming import Session, Field, schema
from ming.declarative import Document
from dbutils import generate_collection_id


__author__ = 'mark'

def _get_dataobj_id():
    """
        Returns the short handle used for this data object in Ming configuration
    """
    return 'simulations'

def _get_collection_id():
    """
    :return: returns the collection name for this data object
    """
    return generate_collection_id("_samples_raw")



def store_stats_mixture_model(config, timestep, num_configs,
                                 config_counts,slatkin,entropy,iqv,
                                 unlabeled_freq, unlabeled_count, conf_slatkin,richness,
                                 unlabeled_freq_ssize,unlabeled_count_ssize,unlabeled_config_ssize,
                                 config_slatkin_ssize,entropy_ssize,iqv_ssize,slatkin_ssize,richness_ssize):
    """Stores the parameters and metadata for a simulation run in the database.


    """
    MixtureModelStats(dict(
        simulation_run_id = config.sim_id,
        sample_time = timestep,
        script_filename = config.script,
        interaction_rule_classes = str(config.INTERACTION_RULE_CLASS),
        pop_class = config.POPULATION_STRUCTURE_CLASS,
        network_class = config.NETWORK_FACTORY_CLASS,
        trait_class = config.TRAIT_FACTORY_CLASS,
        innovation_class = config.INNOVATION_RULE_CLASS,
        num_features = config.num_features,
        init_traits_per_feature = config.num_traits,
        conformism_strength = config.conformism_strength,
        anticonformism_strength = config.anticonformism_strength,
        innovation_rate = config.innovation_rate,
        sample_size = config.SAMPLE_SIZES_STUDIED,
        population_size = config.popsize,
        slatkin_exact = slatkin,
        shannon_entropy = entropy,
        iqv_diversity = iqv,
        num_trait_configurations = num_configs,
        trait_configuration_counts = config_counts,
        unlabeled_frequencies = unlabeled_freq,
        unlabeled_counts = unlabeled_count,
        configuration_slatkin = conf_slatkin,
        pop_richness = richness,
        unlabeled_counts_ssize = unlabeled_count_ssize,
        unlabeled_freq_ssize = unlabeled_freq_ssize,
        unlabeled_config_counts_ssize = unlabeled_config_ssize,
        config_slatkin_ssize = config_slatkin_ssize,
        entropy_ssize = entropy_ssize,
        iqv_ssize = iqv_ssize,
        slatkin_ssize = slatkin_ssize,
        richness_ssize = richness_ssize,
        )).m.insert()
    return True


class MixtureModelStats(Document):

    class __mongometa__:
        session = Session.by_name(_get_dataobj_id())
        name = 'mixture_model_stats'

    _id = Field(schema.ObjectId)
    simulation_run_id = Field(str)
    sample_time = Field(int)
    script_filename = Field(str)
    interaction_rule_classes = Field(str)
    pop_class = Field(str)
    network_class = Field(str)
    trait_class = Field(str)
    innovation_class = Field(str)
    num_features = Field(int)
    init_traits_per_feature = Field(int)
    conformism_strength = Field(float)
    anticonformism_strength = Field(float)
    innovation_rate = Field(float)
    # record all sample sizes used in this simulation run
    sample_size = Field([int])
    population_size = Field(int)
    # whole population statistics
    slatkin_exact = Field([float])
    shannon_entropy = Field([float])
    iqv_diversity = Field([float])
    num_trait_configurations = Field(int)
    trait_configuration_counts = Field([])
    configuration_slatkin = Field(float)
    unlabeled_frequencies = Field([])
    unlabeled_counts = Field([])
    pop_richness = Field([int])
    # results by sample size
    unlabeled_freq_ssize = Field(schema.Anything)
    unlabeled_counts_ssize = Field(schema.Anything)
    unlabeled_config_counts_ssize = Field(schema.Anything)
    num_configurations_ssize = Field(schema.Anything)
    config_slatkin_ssize = Field(schema.Anything)
    entropy_ssize = Field(schema.Anything)
    iqv_ssize = Field(schema.Anything)
    richness_ssize = Field(schema.Anything)
    slatkin_ssize = Field(schema.Anything)




