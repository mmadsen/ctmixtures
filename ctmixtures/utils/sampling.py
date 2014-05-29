#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

import ctmixtures.analysis as analysis
import ctmixtures.data as data


def sample_mixture_model(model, args, config, timestep):
    tfa = analysis.PopulationTraitAnalyzer(model)
    tfa.update()
    ssfa = analysis.SampledTraitAnalyzer(model)
    ssfa.update()


    data.store_stats_mixture_model(config,
                                   timestep,
                                   tfa.get_number_configurations(),
                                   tfa.get_unlabeled_configuration_counts(),
                                   tfa.get_slatkin_exact_probability(),
                                   tfa.get_trait_evenness_entropy(),
                                   tfa.get_trait_evenness_iqv(),
                                   tfa.get_unlabeled_frequency_lists(),
                                   tfa.get_unlableled_count_lists(),
                                   tfa.get_configuration_slatkin_test(),
                                   tfa.get_trait_richness(),
                                   ssfa.get_unlabeled_freq_by_ssize(),
                                   ssfa.get_unlabeled_counts_by_ssize(),
                                   ssfa.get_unlabeled_configuration_counts_by_ssize(),
                                   ssfa.get_configuration_slatkin_by_ssize(),
                                   ssfa.get_entropy_by_ssize(),
                                   ssfa.get_iqv_by_ssize(),
                                   ssfa.get_slatkin_by_ssize(),
                                   ssfa.get_richness_by_ssize()

    )


