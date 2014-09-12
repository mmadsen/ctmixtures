#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Run a single population simulation of a mixture of CT rules, with time averaged observations, trait survival
analysis ala Kandler and Steele, and point/interval measures of trait distribution and diversity.

"""

import logging as log
import argparse
from time import time
import uuid
import numpy.random as npr
import random
import sys

import ming
import ctmixtures.utils as utils
import ctmixtures.data as data
import ctmixtures.analysis as analysis
import pytransmission.popgen as pg


def setup():
    global args, simconfig

    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", help="provide name for experiment", required=True)
    parser.add_argument("--debug", help="turn on debugging output")
    parser.add_argument("--dbhost", help="database hostname, defaults to localhost", default="localhost")
    parser.add_argument("--dbport", help="database port, defaults to 27017", default="27017")
    parser.add_argument("--configuration", help="Configuration file for experiment", required=True)
    parser.add_argument("--popsize", help="Population size", required=True)
    parser.add_argument("--numloci", help="Number of loci per individual", required=True)
    parser.add_argument("--maxinittraits", help="Max initial number of traits per locus for initialization", required=True)
    parser.add_argument("--conformismstrength", help="Strength of conformist bias [0.0 - 1.0]")
    parser.add_argument("--anticonformismstrength", help="Strength of conformist bias [0.0 - 1.0]")
    parser.add_argument("--innovationrate", help="Theta value rate at which innovations occur in population", required=True)
    parser.add_argument("--periodic", help="Periodic boundary condition", choices=['1','0'], required=True)
    parser.add_argument("--kandlerinterval", help="Interval for Kandler remaining traits sample, taken before maxtime, in generations (will be scaled to timesteps)", default="1000")
    parser.add_argument("--simulationendtime", help="Time at which simulation and sampling end, defaults to 2M steps", default="2000000")
    parser.add_argument("--seed", type=int, help="Seed for random generators to ensure replicability")

    args = parser.parse_args()

    simconfig = utils.MixtureConfiguration(args.configuration)

    if args.debug == '1':
        log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    else:
        log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')

    if args.seed is None:
        log.debug("No seed given, allowing RNG's to initialize randomly")
    else:
        log.debug("Seeding RNGs with seed: %s", args.seed)
        npr.seed(args.seed)
        random.seed(args.seed)
        simconfig.random_seed = args.seed

    simconfig.full_command_line = " ".join(sys.argv)

    log.debug("experiment name: %s", args.experiment)
    data.set_experiment_name(args.experiment)
    data.set_database_hostname(args.dbhost)
    data.set_database_port(args.dbport)
    config = data.getMingConfiguration(data.modules)
    ming.configure(**config)

    simconfig.num_features = int(args.numloci)
    simconfig.num_traits = int(args.maxinittraits)
    simconfig.popsize = int(args.popsize)
    mut = pg.moran_mutation_rate_from_theta(float(args.popsize), float(args.innovationrate))
    simconfig.innovation_rate = float(args.numloci) * mut
    simconfig.configured_innovation_rate = float(args.innovationrate)
    log.debug("configured theta = %s, using numloci %s * per-locus mutation rate %s = all-loci innovation rate: %s", args.innovationrate, args.numloci, mut, simconfig.innovation_rate)


    simconfig.maxtime = int(args.simulationendtime)
    simconfig.script = __file__
    simconfig.conformism_strength = float(args.conformismstrength)
    simconfig.anticonformism_strength = float(args.anticonformismstrength)
    simconfig.maxtime = int(args.simulationendtime)
    simconfig.model_class_label = simconfig.MODEL_CLASS_LABEL

    log.debug("Equifinality model class: %s", simconfig.model_class_label)


    simconfig.sim_id = uuid.uuid4().urn
    if args.periodic == '1':
        simconfig.periodic = 1
    else:
        simconfig.periodic = 0


def main():
    start = time()
    log.debug("Configuring CT Mixture Model with structure class: %s graph factory: %s interaction rule: %s", simconfig.POPULATION_STRUCTURE_CLASS, simconfig.NETWORK_FACTORY_CLASS, simconfig.INTERACTION_RULE_CLASS)


    # calculate timing and intervals for samples, given parameters
    ta_interval_list = simconfig.TIME_AVERAGING_DURATIONS
    max_ta_interval = max(ta_interval_list)
    kandler_interval_in_generations = int(args.kandlerinterval)
    kandler_interval_timesteps = kandler_interval_in_generations * simconfig.popsize
    ending_indextime = simconfig.maxtime - (max_ta_interval * simconfig.popsize)
    starting_indextime = ending_indextime - kandler_interval_timesteps

    log.debug("max TA interval: %s", max_ta_interval)
    log.debug("Kandler survival interval in generations: %s  timesteps: %s", kandler_interval_in_generations, kandler_interval_timesteps)
    log.debug("starting Kandler survival interval tick: %s", starting_indextime)
    log.debug("starting Kandler survival interval tick: %s", ending_indextime)

    perlocus_innovation_rate = pg.moran_mutation_rate_from_theta(simconfig.popsize, simconfig.configured_innovation_rate)
    multilocus_stationarity_time = pg.moran_watkins_multilocus_convergence_time_timesteps(simconfig.popsize, simconfig.num_features, perlocus_innovation_rate)

    model_constructor = utils.load_class(simconfig.POPULATION_STRUCTURE_CLASS)
    graph_factory_constructor = utils.load_class(simconfig.NETWORK_FACTORY_CLASS)
    trait_factory_constructor = utils.load_class(simconfig.TRAIT_FACTORY_CLASS)
    interaction_rule_list = utils.parse_interaction_rule_map(simconfig.INTERACTION_RULE_CLASS)
    innovation_rule_constructor = utils.load_class(simconfig.INNOVATION_RULE_CLASS)
    dynamics_constructor = utils.load_class(simconfig.DYNAMICS_CLASS)
    timeaveraging_constructor = utils.load_class(simconfig.TIME_AVERAGING_CLASS)

    # instantiate the model and its various subobjects, including any interaction rules
    graph_factory = graph_factory_constructor(simconfig)
    trait_factory = trait_factory_constructor(simconfig)


    # construct two time averaging objects
    starting_ta_sample = timeaveraging_constructor(starting_indextime, ta_interval_list, simconfig.popsize, simconfig.num_features, ending_interval=False)
    ending_ta_sample = timeaveraging_constructor(ending_indextime, ta_interval_list, simconfig.popsize, simconfig.num_features, ending_interval=True)

    # calculate times for various sampling events
    earliest_sample_time = starting_ta_sample.get_earliest_tick_for_all_intervals()
    kandler_start_time_nota = starting_ta_sample.get_latest_tick_for_all_intervals()
    kandler_stop_time_nota = ending_ta_sample.get_earliest_tick_for_all_intervals()

    log.debug("starting TA sampler intervals: %s", starting_ta_sample.get_interval_tuples())
    log.debug("ending TA sampler intervals: %s", ending_ta_sample.get_interval_tuples())
    log.debug("earliest sample time: %s", earliest_sample_time)
    log.debug("kandler interval start time: %s", kandler_start_time_nota)
    log.debug("kandler interval stop time: %s", kandler_stop_time_nota)

    log.debug("Minimum stationarity time: %s", multilocus_stationarity_time)

    # if multilocus_stationarity_time > earliest_sample_time:
    #     log.error("Sampling start time is before stationarity reached, increase simulation max time!")
    #     exit(0)


    model = model_constructor(simconfig, graph_factory, trait_factory)
    interaction_rule_list = utils.construct_rule_objects(interaction_rule_list, model)
    model.interaction_rules = interaction_rule_list

    # now we're ready to initialize the population
    model.initialize_population()
    innovation_rule = innovation_rule_constructor(model)

    # initialize a dynamics
    dynamics = dynamics_constructor(simconfig,model,innovation_rule)

    tfa = analysis.TimeAveragedSampledTraitAnalyzer(model, starting_ta_sample, ending_ta_sample)
    ssfa = analysis.SampledTraitAnalyzer(model)


    log.info("Kandler tracking interval start: %s  stop: %s", kandler_start_time_nota, kandler_stop_time_nota)


    log.info("Starting %s", simconfig.sim_id)

    while(1):

        timestep = dynamics.update()

        if (timestep % 10000) == 0:
            log.debug("time: %s  copies by locus: %s  innovations: %s innov by locus: %s",
                      timestep, model.get_interactions_by_locus(), model.get_innovations(),
                      model.get_innovations_by_locus())

        # starting at the beginning of the first time averaging window, we start feeding updates
        # to the time averagers

        if (timestep == earliest_sample_time):
            log.info("Hit earliest sample time:  %s", timestep)

        if ( timestep >= earliest_sample_time ):
            tfa.update(timestep)

        # We want to sample the whole population for a Kandler & Shennan style trait survival analysis.
        # This is the full population synchronic calculation, so it happens at two points in time.
        # otherwise, the time averagers take care of the case where start & end are non-trivial duration samples.
        if ( timestep == kandler_start_time_nota ):
            utils.start_kandler_remaining_trait_tracking(tfa, ssfa, timestep)

        if ( timestep == kandler_stop_time_nota ):
            utils.stop_kandler_remaining_trait_tracking(tfa, ssfa, timestep)


        # sample and end the simulation, also recording time averaged stats for the population,
        # and kandler & shennan trait survival for the time averaged samples.
        # finally, record simulation timing so we can track performance and plan blocks of simulation runs
        if timestep >= simconfig.maxtime:
            utils.record_final_samples(tfa, ssfa, simconfig, timestep)
            endtime = time()
            elapsed = endtime - start
            log.info("Completed: %s  Elapsed: %s", simconfig.sim_id, elapsed)
            data.store_simulation_timing(simconfig.sim_id,simconfig.INTERACTION_RULE_CLASS,simconfig.POPULATION_STRUCTURE_CLASS,simconfig.script,args.experiment,elapsed,timestep,simconfig.popsize)
            exit(0)

# end main




if __name__ == "__main__":
    setup()
    main()

