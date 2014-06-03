#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
the purpose here is to run the tree structured model for a fixed length of time, with mutation, and see what
happens to the structure of trait trees, sampled at fixed intervals.  We will leave in homophily, but ignore
convergence since it should be temporary with mutation/noise.

"""

import logging as log
import argparse
from time import time
import uuid

import ming

import ctmixtures.utils as utils
import ctmixtures.data as data
import ctmixtures.dynamics as dyn
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
    parser.add_argument("--conformismstrength", help="Strength of conformist bias [0.0 - 1.0]", required=True)
    parser.add_argument("--anticonformismstrength", help="Strength of conformist bias [0.0 - 1.0]", required=True)
    parser.add_argument("--innovationrate", help="Theta value rate at which innovations occur in population", required=True)
    parser.add_argument("--periodic", help="Periodic boundary condition", choices=['1','0'], required=True)
    parser.add_argument("--kandlerinterval", help="Interval for Kandler remaining traits sample, taken before maxtime, in generations (will be scaled to timesteps)", default="1000")
    parser.add_argument("--simulationendtime", help="Time at which simulation and sampling end, defaults to 2M steps", default="2000000")

    args = parser.parse_args()

    simconfig = utils.MixtureConfiguration(args.configuration)

    if args.debug == '1':
        log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    else:
        log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')

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
    log.debug("configured theta = %s, using numloci %s * per-locus mutation rate %s = all-loci innovation rate: %s", args.innovationrate, args.numloci, mut, simconfig.innovation_rate)


    simconfig.maxtime = int(args.simulationendtime)
    simconfig.script = __file__
    simconfig.conformism_strength = float(args.conformismstrength)
    simconfig.anticonformism_strength = float(args.anticonformismstrength)
    simconfig.maxtime = int(args.simulationendtime)


    simconfig.sim_id = uuid.uuid4().urn
    if args.periodic == '1':
        simconfig.periodic = 1
    else:
        simconfig.periodic = 0


def main():
    start = time()

    kandler_interval_in_generations = int(args.kandlerinterval)
    kandler_interval_timesteps = kandler_interval_in_generations * simconfig.popsize
    kandler_start_time = simconfig.maxtime - kandler_interval_timesteps

    log.debug("Taking a Kandler trait survival sample of %s timesteps, beginning at tick %s", kandler_interval_timesteps, kandler_start_time)

    model_constructor = utils.load_class(simconfig.POPULATION_STRUCTURE_CLASS)
    graph_factory_constructor = utils.load_class(simconfig.NETWORK_FACTORY_CLASS)
    trait_factory_constructor = utils.load_class(simconfig.TRAIT_FACTORY_CLASS)
    interaction_rule_list = utils.parse_interaction_rule_map(simconfig.INTERACTION_RULE_CLASS)
    innovation_rule_constructor = utils.load_class(simconfig.INNOVATION_RULE_CLASS)

    log.debug("Configuring CT Mixture Model with structure class: %s graph factory: %s interaction rule: %s", simconfig.POPULATION_STRUCTURE_CLASS, simconfig.NETWORK_FACTORY_CLASS, simconfig.INTERACTION_RULE_CLASS)

    # instantiate the model and its various subobjects, including any interaction rules
    graph_factory = graph_factory_constructor(simconfig)
    trait_factory = trait_factory_constructor(simconfig)



    model = model_constructor(simconfig, graph_factory, trait_factory)
    interaction_rule_list = utils.construct_rule_objects(interaction_rule_list, model)
    model.interaction_rules = interaction_rule_list

    # now we're ready to initialize the population
    model.initialize_population()
    innovation_rule = innovation_rule_constructor(model)

    # initialize a dynamics
    dynamics = dyn.MoranDynamics(simconfig,model,innovation_rule)


    tfa = analysis.PopulationTraitAnalyzer(model)
    ssfa = analysis.SampledTraitAnalyzer(model)

    log.info("Starting %s", simconfig.sim_id)

    #if (args.debug == '1'):
        #utils.debug_sample_mixture_model(tfa, ssfa, simconfig, 0)

    while(1):

        dynamics.update()
        timestep = dynamics.timestep

        if (timestep % 100000) == 0:
            log.debug("time: %s  copies by locus: %s  innovations: %s innov by locus: %s",
                      timestep, model.get_interactions_by_locus(), model.get_innovations(),
                      model.get_innovations_by_locus())
            #utils.debug_sample_mixture_model(tfa, ssfa, simconfig, timestep)

        if ( timestep == kandler_start_time ):
            utils.start_kandler_remaining_trait_tracking(tfa, ssfa, timestep)

        # sample and end the simulation
        if timestep >= simconfig.maxtime:
            utils.stop_kandler_remaining_trait_tracking(tfa, ssfa, timestep)
            utils.sample_mixture_and_record_kandler_remaining_traits(tfa, ssfa, simconfig, timestep)
            endtime = time()
            elapsed = endtime - start
            log.info("Completed: %s  Elapsed: %s", simconfig.sim_id, elapsed)
            data.store_simulation_timing(simconfig.sim_id,simconfig.INTERACTION_RULE_CLASS,simconfig.POPULATION_STRUCTURE_CLASS,simconfig.script,args.experiment,elapsed,timestep)
            exit(0)

# end main




if __name__ == "__main__":
    setup()
    main()

