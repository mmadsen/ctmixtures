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
import numpy.random as npr
import json


def generate_random_neutral_model(seed):
    """
    Creates a simulation command line for the sim-ctmixture-timeaveraging.py model, using a random
    value for --innovationrate drawn uniformly from a configured prior range.

    :return: string
    """

    theta = npr.uniform(low = float(expconfig['theta_low']), high = float(expconfig['theta_high']))

    if args.simprefix is not None:
        cmd = ""
        cmd += args.simprefix
        cmd += "/sim-ctmixture-timeaveraging.py "
    else:
        cmd = "sim-ctmixture-timeaveraging.py "

    cmd += " --experiment "
    cmd += args.experiment
    cmd += " --configuration "

    if args.confprefix is not None:
        cmd += args.confprefix
        cmd += "/"
        cmd += args.configuration
    else:
        cmd += args.configuration

    cmd += " --popsize "
    cmd += str(expconfig["popsize"])
    cmd += " --maxinittraits "
    cmd += str(expconfig["maxinittraits"])
    cmd += " --numloci "
    cmd += str(expconfig["numloci"])
    cmd += " --kandlerinterval "
    cmd += str(expconfig["kandler_interval"])
    cmd += " --innovationrate "
    cmd += str(theta)
    cmd += " --periodic 0 "
    cmd += " --simulationendtime "
    cmd += str(expconfig["endtime"])
    cmd += " --conformismstrength "
    cmd += str(0.0)
    cmd += " --anticonformismstrength "
    cmd += str(0.0)
    cmd += " --debug "
    cmd += args.debug
    cmd += " --seed "
    cmd += str(seed)
    cmd += '\n'

    #log.debug("%s", cmd)
    return cmd


def generate_random_conformist_model(seed):
    """
    Creates a simulation command line for the sim-ctmixture-timeaveraging.py model, using a random
    value for --innovationrate drawn uniformly from a configured prior range.

    :return: string
    """

    theta = npr.uniform(low = float(expconfig['theta_low']), high = float(expconfig['theta_high']))
    conf_str = npr.uniform(low = float(expconfig['conformist_strength_low']), high = float(expconfig['conformist_strength_high']))
    aconf_str = npr.uniform(low = float(expconfig['anticonformist_strength_low']), high = float(expconfig['anticonformist_strength_high']))

    if args.simprefix is not None:
        cmd = ""
        cmd += args.simprefix
        cmd += "/sim-ctmixture-timeaveraging.py "
    else:
        cmd = "sim-ctmixture-timeaveraging.py "

    cmd += " --experiment "
    cmd += args.experiment
    cmd += " --configuration "

    if args.confprefix is not None:
        cmd += args.confprefix
        cmd += "/"
        cmd += args.configuration
    else:
        cmd += args.configuration

    cmd += " --popsize "
    cmd += str(expconfig["popsize"])
    cmd += " --maxinittraits "
    cmd += str(expconfig["maxinittraits"])
    cmd += " --numloci "
    cmd += str(expconfig["numloci"])
    cmd += " --kandlerinterval "
    cmd += str(expconfig["kandler_interval"])
    cmd += " --innovationrate "
    cmd += str(theta)
    cmd += " --periodic 0 "
    cmd += " --simulationendtime "
    cmd += str(expconfig["endtime"])
    cmd += " --conformismstrength "
    cmd += str(conf_str)
    cmd += " --anticonformismstrength "
    cmd += str(aconf_str)
    cmd += " --debug "
    cmd += args.debug
    cmd += " --seed "
    cmd += str(seed)

    cmd += '\n'

    #log.debug("%s", cmd)
    return cmd


def parse_experiment_config(config_path):
    try:
        json_data = open(config_path)
        config = json.load(json_data)
    except ValueError:
        print "Problem parsing json configuration file - probably malformed syntax"
        exit(1)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
        exit(1)

    return config



def setup():
    global args, simconfig, expconfig

    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", help="provide name for experiment", required=True)
    parser.add_argument("--expconfig", help="Experiment configuration file path", required=True)
    parser.add_argument("--debug", help="turn on debugging output")
    parser.add_argument("--dbhost", help="database hostname, defaults to localhost", default="localhost")
    parser.add_argument("--dbport", help="database port, defaults to 27017", default="27017")
    parser.add_argument("--configuration", help="Configuration file to use in generated simulation commands", required=True)
    parser.add_argument("--parallelism", help="Number of separate generated command lists to create", default="4")
    parser.add_argument("--numsims", type=int, help="Number of simulations to generate by random prior sampling", default=100)
    parser.add_argument("--model", help="Model being generated", choices=['neutral', 'conformist'])
    parser.add_argument("--simprefix", help="Full path prefix to the simulation executable (optional)")
    parser.add_argument("--confprefix", help="Full path prefix to the configuration file (optional)")

    args = parser.parse_args()

    simconfig = utils.MixtureConfiguration(args.configuration)
    expconfig = parse_experiment_config(args.expconfig)

    if args.debug == '1':
        log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    elif args.debug is None:
        args.debug = '0'
    else:
        log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')

    log.info("Generating simulation commands for experiment: %s and model: %s", args.experiment, args.model)



def main():
    log.info("Opening %s output files for simulation configuration", args.parallelism)
    num_files = int(args.parallelism)
    file_list = []
    base_name = "simrunner-exp-"
    base_name += args.experiment
    base_name += "-"
    base_name += args.model
    base_name += "-"

    for i in range(0, num_files):
        filename = ''
        filename += base_name
        filename += str(i)
        filename += ".sh"

        f = open(filename, 'w')

        f.write("#!/bin/sh\n\n")
        file_list.append(f)

    file_cycle = itertools.cycle(file_list)


    for i in xrange(0, args.numsims):

        # give us a random seed that will fit in a 64 bit long integer
        seed = npr.randint(1,2**62)

        if args.model == 'neutral':
            cmd = generate_random_neutral_model(seed)
        elif args.model == 'conformist':
            cmd = generate_random_conformist_model(seed)
        else:
            log.error("unrecognized model type, add to argparse and create a generator method")
            exit(1)


        fc = file_cycle.next()
        fc.write(cmd)


    for fh in file_list:
        fh.close()





if __name__ == "__main__":
    setup()
    main()


