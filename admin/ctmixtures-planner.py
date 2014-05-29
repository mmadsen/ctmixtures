#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

import logging as log
import argparse
import itertools

import ctmixtures.utils as utils


def setup():
    global args, simconfig

    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", help="provide name for experiment", required=False)
    parser.add_argument("--debug", help="turn on debugging output")
    parser.add_argument("--dbhost", help="database hostname, defaults to localhost", default="localhost")
    parser.add_argument("--dbport", help="database port, defaults to 27017", default="27017")
    parser.add_argument("--configuration", help="Configuration file for experiment", required=True)


    args = parser.parse_args()

    simconfig = utils.MixtureConfiguration(args.configuration)

    if args.debug == '1':
        log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    else:
        log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')

    log.debug("experiment name: %s", args.experiment)



def main():

    structure_class_name = simconfig.POPULATION_STRUCTURE_CLASS
    log.info("Configuring CTMixture model with structure class: %s", structure_class_name)


    basic_config = utils.MixtureConfiguration(args.configuration)

    state_space = [
        basic_config.POPULATION_SIZES_STUDIED,
        basic_config.INNOVATION_RATE,
        basic_config.NUMBER_OF_DIMENSIONS_OR_FEATURES,
        basic_config.NUMBER_OF_TRAITS_PER_DIMENSION,
        basic_config.CONFORMISM_STRENGTH,
        basic_config.ANTICONFORMISM_STRENGTH
    ]


    if basic_config.NETWORK_FACTORY_CLASS == 'madsenlab.ctmixtures.population.WattsStrogatzSmallWorldFactory':
        state_space.append(basic_config.WS_REWIRING_FACTOR)



    num_runs = 0

    for param_combination in itertools.product(*state_space):
        for replication in range(0, basic_config.REPLICATIONS_PER_PARAM_SET):
            num_runs += 1


    log.info("Total number of runs: %s", num_runs)


if __name__ == "__main__":
    setup()
    main()


