#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

import ctmixtures.analysis as analysis
import ctmixtures.data as data
import logging as log



def debug_sample_mixture_model(tfa, ssfa, config, timestep):
    tfa.update(timestep)
    ssfa.update(timestep)
    log.debug("trait counts: %s", tfa.get_trait_counts())


def sample_mixture_model(tfa, ssfa, config, timestep):
    tfa.update(timestep)
    ssfa.update(timestep)


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
                                   ssfa.get_richness_by_ssize(),
                                   None, #kandler interval
                                   None, #kandler_remaining_count
                                   None, # TAssize fields follow
                                   None,
                                   None,
                                   None,
                                   None,
                                   None,
                                   None,
                                   None,
                                   None,
                                   None,
                                   None
    )



def start_kandler_remaining_trait_tracking(tfa, ssfa, timestep):
    log.debug("Starting Kandler remaining trait tracking at %s", timestep)
    tfa.update(timestep)
    tfa.kandler_survival_start(timestep)

def stop_kandler_remaining_trait_tracking(tfa, ssfa, timestep):
    log.debug("Stopping Kandler remaining trait tracking at %s", timestep)
    tfa.update(timestep)
    tfa.kandler_survival_stop(timestep)

def record_final_samples(tfa, ssfa, config, timestep):
    tfa.update(timestep)
    ssfa.update(timestep)
    # sample from the time averagers, for later statistics calculation
    # only applies to the TA and sampled analyzer, the others will simply log an
    # unknown method call
    tfa.take_sample_snapshot()

    (interval, remaining_traits) = tfa.get_kandler_remaining_traits_per_locus()


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
                                   ssfa.get_richness_by_ssize(),
                                   interval,
                                   remaining_traits,
                                   tfa.get_ta_unlabeled_frequency_lists(),
                                   tfa.get_ta_trait_richness(),
                                   tfa.get_ta_slatkin_exact_probability(),
                                   tfa.get_ta_trait_evenness_entropy(),
                                   tfa.get_ta_trait_evenness_iqv(),
                                   tfa.get_ta_unlabeled_configuration_counts(),
                                   tfa.get_ta_number_configurations(),
                                   tfa.get_ta_configuration_slatkin_test(),
                                   tfa.get_ta_configuration_evenness_entropy(),
                                   tfa.get_ta_configuration_evenness_iqv(),
                                   tfa.get_ta_kandler_remaining_traits_per_locus()
    )


