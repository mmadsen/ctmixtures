#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""
import logging as log

from configuration import MixtureConfiguration
from dynamicloading import load_class, parse_interaction_rule_map, construct_rule_objects
from ctmixtures.utils.sampling import sample_mixture_model
import pytransmission.popgen as pg


def calc_simulation_endtime(config, popsize, innovation_rate, num_loci):
    """
    We want to end the simulation at a reasonable point after we have reached
    quasi-stationary equilibrium and taken several samples.

    :param popsize:
    :param innovation_rate:
    :param num_loci:
    :return:
    """
    stationarity_in_ticks = pg.moran_watkins_multilocus_convergence_time_timesteps(popsize,
                                                                                   innovation_rate,
                                                                                   num_loci)

    log.debug("calc_sim_endtime - stationarity: %s", stationarity_in_ticks)

