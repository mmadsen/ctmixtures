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

from ming import Session, Field, schema
from ming.declarative import Document

from ctmixtures.data.dbutils import generate_collection_id


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
                                 config_slatkin_ssize,entropy_ssize,iqv_ssize,slatkin_ssize,richness_ssize,
                                 kandler_interval,kandler_remaining_count,unlab_freq_tassize,
                                 richness_tassize,slatkin_tassize,entropy_tassize,iqv_tassize,
                                 unlab_ccount_tassize,config_richness_tassize,config_slatkin_tassize,
                                 config_entropy_tassize,config_iqv_tassize,kandler_remaining_tassize):
    """Stores the parameters and metadata for a simulation run in the database.
    """
    MixtureModelStats(dict(
        simulation_run_id = config.sim_id,
        sample_time = timestep,
        script_filename = config.script,
        model_class_label = config.model_class_label,
        interaction_rule_classes = str(config.INTERACTION_RULE_CLASS),
        pop_class = config.POPULATION_STRUCTURE_CLASS,
        network_class = config.NETWORK_FACTORY_CLASS,
        trait_class = config.TRAIT_FACTORY_CLASS,
        innovation_class = config.INNOVATION_RULE_CLASS,
        num_features = config.num_features,
        init_traits_per_feature = config.num_traits,
        conformism_strength = config.conformism_strength,
        anticonformism_strength = config.anticonformism_strength,
        innovation_rate = config.configured_innovation_rate,
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
        kandler_interval = kandler_interval,
        kandler_remaining_count = kandler_remaining_count,
        unlabeled_freq_ta_ssize = unlab_freq_tassize,
        richness_ta_ssize = richness_tassize,
        slatkin_ta_ssize = slatkin_tassize,
        entropy_ta_ssize = entropy_tassize,
        iqv_ta_ssize = iqv_tassize,
        unlabeled_config_counts_ta_ssize = unlab_ccount_tassize,
        num_configurations_ta_ssize = config_richness_tassize,
        config_slatkin_ta_ssize = config_slatkin_tassize,
        config_entropy_ta_ssize = config_entropy_tassize,
        config_iqv_ta_ssize = config_iqv_tassize,
        kandler_remaining_tassize = kandler_remaining_tassize
        )).m.insert()
    return True


def columns_to_export_for_analysis():
    """

    :return:
    """
    cols = [
        "simulation_run_id",
        "model_class_label",
        "network_class",
        "interaction_rule_classes",
        "sample_time",
        "num_features",
        "population_size",
        "innovation_rate",
        "conformism_strength",
        "anticonformism_strength",
        "kandler_interval"
    ]
    return cols

def tassize_columns_to_export_for_analysis():
    cols = [
        "simulation_run_id",
        "network_class",
        "interaction_rule_classes",
        "sample_time",
        "num_features",
        "population_size",
        "innovation_rate",
        "conformism_strength",
        "anticonformism_strength",
        "kandler_interval"
    ]

def ssize_columns_to_export_for_analysis():
    cols = [
        "simulation_run_id",
        "network_class",
        "interaction_rule_classes",
        "sample_time",
        "num_features",
        "population_size",
        "innovation_rate",
        "conformism_strength",
        "anticonformism_strength",
        "kandler_interval"
    ]





class MixtureModelStats(Document):

    class __mongometa__:
        session = Session.by_name(_get_dataobj_id())
        name = 'mixture_model_stats'

    # model and parameters
    _id = Field(schema.ObjectId)
    simulation_run_id = Field(str)
    sample_time = Field(int)
    model_class_label = Field(str)
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
    sample_size = Field([int])
    population_size = Field(int)
    kandler_interval = Field(int)

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
    kandler_remaining_count = Field([int])


    # results for TA intervals over all sample sizes
    unlabeled_freq_ta_ssize = Field(schema.Anything)
    richness_ta_ssize = Field(schema.Anything)
    slatkin_ta_ssize = Field(schema.Anything)
    entropy_ta_ssize = Field(schema.Anything)
    iqv_ta_ssize = Field(schema.Anything)
    unlabeled_config_counts_ta_ssize = Field(schema.Anything)
    num_configurations_ta_ssize = Field(schema.Anything)
    config_slatkin_ta_ssize = Field(schema.Anything)
    config_entropy_ta_ssize = Field(schema.Anything)
    config_iqv_ta_ssize = Field(schema.Anything)
    kandler_remaining_tassize = Field(schema.Anything)


# TODO - add final set of fields to storage function above


