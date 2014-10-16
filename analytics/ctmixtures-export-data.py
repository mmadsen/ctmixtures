#!/usr/bin/env python

# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

import ming
import csv
import logging as log
import argparse
import ctmixtures.data as data
import numpy as np


############################################################################
def setup():
    global args, config, simconfig
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", help="provide name for experiment, to be used as prefix for database collections")
    parser.add_argument("--debug", help="turn on debugging output")
    parser.add_argument("--dbhost", help="database hostname, defaults to localhost", default="localhost")
    parser.add_argument("--dbport", help="database port, defaults to 27017", default="27017")
    parser.add_argument("--configuration", help="Path to configuration file")
    parser.add_argument("--filename", help="path and base filename for exports (DO NOT include *.csv extension)", required=True)

    args = parser.parse_args()

    if args.debug == 1:
        log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    else:
        log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')

    #### main program ####
    log.info("EXPORT DATA TO CSV - Experiment: %s", args.experiment)
    data.set_experiment_name(args.experiment)
    data.set_database_hostname(args.dbhost)
    data.set_database_port(args.dbport)
    config = data.getMingConfiguration(data.modules)
    ming.configure(**config)



############################################################################
def export_simulation_record():
    # ## Export a simulation record file, with all params and classes used, random
    ### seed, whatever is needed to replicate the simulations
    full_filename = ''
    full_filename += args.filename
    full_filename += "-simulation-data.csv"
    sim_fields = data.mixture_model_stats.sim_record_columns_to_export()
    ofile = open(full_filename, "wb")
    writer = csv.DictWriter(ofile, fieldnames=sim_fields, quotechar='"', quoting=csv.QUOTE_ALL)
    headers = dict((n, n) for n in sim_fields)
    writer.writerow(headers)
    cursor = data.MixtureModelStats.m.find(dict(), dict(timeout=False))
    for sample in cursor:
        row = dict()
        for field in sim_fields:
            row[field] = sample[field]

        # correct kandler_interval from timesteps to generations
        row['kandler_interval'] = int(row['kandler_interval']) / int(row['population_size'])

        #log.info("sim data row: %s", row)
        writer.writerow(row)
    ofile.close()


############################################################################
# # whole population statistics
# slatkin_exact = Field([float])
# shannon_entropy = Field([float])
# iqv_diversity = Field([float])
# num_trait_configurations = Field(int)
# trait_configuration_counts = Field([])
# configuration_slatkin = Field(float)
# unlabeled_frequencies = Field([])
# unlabeled_counts = Field([])
# pop_richness = Field([int])

def export_population_stats():
    # ## Export a full population census statistics file ###
    full_filename = ''
    full_filename += args.filename
    full_filename += "-population-data.csv"
    pop_fields = data.mixture_model_stats.pop_columns_to_export()

    # adjust the fields for the new summary statistics
    pop_fields.append('innovation_rate')
    pop_fields.append('slatkin_locus_max')
    pop_fields.append('slatkin_locus_min')
    pop_fields.append('entropy_locus_max')
    pop_fields.append('entropy_locus_min')
    pop_fields.append('entropy_locus_mean')
    pop_fields.append('iqv_locus_max')
    pop_fields.append('iqv_locus_min')
    pop_fields.append('iqv_locus_mean')
    pop_fields.append('richness_locus_max')
    pop_fields.append('richness_locus_min')
    pop_fields.append('richness_locus_mean')
    pop_fields.append('kandler_locus_max')
    pop_fields.append('kandler_locus_min')
    pop_fields.append('kandler_locus_mean')

    ofile = open(full_filename, "wb")
    writer = csv.DictWriter(ofile, fieldnames=pop_fields, quotechar='"', quoting=csv.QUOTE_ALL)
    headers = dict((n, n) for n in pop_fields)
    writer.writerow(headers)

    cursor = data.MixtureModelStats.m.find(dict(), dict(timeout=False))
    for sample in cursor:
        row = dict()
        row['simulation_run_id'] = sample['simulation_run_id']
        row['model_class_label'] = sample['model_class_label']
        row['num_trait_configurations'] = sample['num_trait_configurations']
        row['configuration_slatkin'] = sample['configuration_slatkin']
        row['innovation_rate'] = sample['innovation_rate']

        # slatkin exact
        slatkin_values = sample['slatkin_exact']
        row['slatkin_locus_max'] = max(slatkin_values)
        row['slatkin_locus_min'] = min(slatkin_values)

        # shannon entropy
        entropy_list = sample['slatkin_exact']
        row['entropy_locus_max'] = max(entropy_list)
        row['entropy_locus_min'] = min(entropy_list)
        row['entropy_locus_mean'] = np.average(entropy_list)

        # IQV
        iqv_list = sample['iqv_diversity']
        row['iqv_locus_max'] = max(iqv_list)
        row['iqv_locus_min'] = min(iqv_list)
        row['iqv_locus_mean'] = np.average(iqv_list)

        # Per-locus richness
        richness_list = sample['pop_richness']
        row['richness_locus_max'] = max(richness_list)
        row['richness_locus_min'] = min(richness_list)
        row['richness_locus_mean'] = np.average(richness_list)

        # Kandler remaining per locus
        kandler_list = sample['kandler_remaining_count']
        row['kandler_locus_max'] = max(kandler_list)
        row['kandler_locus_min'] = min(kandler_list)
        row['kandler_locus_mean'] = np.average(kandler_list)

        #log.info("sim data row: %s", row)
        writer.writerow(row)
    ofile.close()

############################################################################
# # results by sample size
# unlabeled_freq_ssize = Field(schema.Anything)
# unlabeled_counts_ssize = Field(schema.Anything)
# unlabeled_config_counts_ssize = Field(schema.Anything)
# num_configurations_ssize = Field(schema.Anything)
# config_slatkin_ssize = Field(schema.Anything)
# entropy_ssize = Field(schema.Anything)
# iqv_ssize = Field(schema.Anything)
# richness_ssize = Field(schema.Anything)
# slatkin_ssize = Field(schema.Anything)
# kandler_remaining_count = Field([int])

def export_sampled_stats():
    ## export a file with sampled statistics
    full_filename = ''
    full_filename += args.filename
    full_filename += "-sampled-data.csv"

    sim_fields = data.mixture_model_stats.ssize_columns_to_export()

    sim_fields.append('innovation_rate')
    sim_fields.append('sample_size')
    sim_fields.append('config_slatkin_ssize')
    sim_fields.append('num_configurations_ssize')
    sim_fields.append('richness_locus_min')
    sim_fields.append('richness_locus_max')
    sim_fields.append('richness_locus_mean')
    sim_fields.append('slatkin_locus_min')
    sim_fields.append('slatkin_locus_max')
    sim_fields.append('slatkin_locus_mean')
    sim_fields.append('entropy_locus_min')
    sim_fields.append('entropy_locus_max')
    sim_fields.append('entropy_locus_mean')
    sim_fields.append('iqv_locus_min')
    sim_fields.append('iqv_locus_max')
    sim_fields.append('iqv_locus_mean')

    ofile = open(full_filename, "wb")
    writer = csv.DictWriter(ofile, fieldnames=sim_fields, quotechar='"', quoting=csv.QUOTE_ALL)
    headers = dict((n, n) for n in sim_fields)
    writer.writerow(headers)
    cursor = data.MixtureModelStats.m.find(dict(), dict(timeout=False))

    for sample in cursor:
        row = dict()
        row['simulation_run_id'] = sample['simulation_run_id']
        row['model_class_label'] = sample['model_class_label']
        row['innovation_rate'] = sample['innovation_rate']

        # all of the other fields require iterating over sample sizes.
        ssizes = sample['sample_size']
        for ssize in ssizes:
            row['sample_size'] = ssize
            row['config_slatkin_ssize'] = sample['config_slatkin_ssize'][str(ssize)]

            config_cnt = sample['unlabeled_config_counts_ssize'][str(ssize)]
            row['num_configurations_ssize'] = len(config_cnt)

            richness_list = sample['richness_ssize'][str(ssize)]
            row['richness_locus_min'] = min(richness_list)
            row['richness_locus_max'] = max(richness_list)
            row['richness_locus_mean'] = np.average(richness_list)

            slatkin_list = sample['slatkin_ssize'][str(ssize)]
            row['slatkin_locus_min'] = min(slatkin_list)
            row['slatkin_locus_max'] = max(slatkin_list)
            row['slatkin_locus_mean'] = np.average(slatkin_list)

            entropy_list = sample['entropy_ssize'][str(ssize)]
            row['entropy_locus_min'] = min(entropy_list)
            row['entropy_locus_max'] = max(entropy_list)
            row['entropy_locus_mean'] = np.average(entropy_list)

            iqv_list = sample['iqv_ssize'][str(ssize)]
            row['iqv_locus_min'] = min(iqv_list)
            row['iqv_locus_max'] = max(iqv_list)
            row['iqv_locus_mean'] = np.average(iqv_list)

        #log.info("sim data row: %s", row)
        writer.writerow(row)

    ofile.close()

############################################################################
# # results for TA intervals over all sample sizes
# unlabeled_freq_ta_ssize = Field(schema.Anything)
# richness_ta_ssize = Field(schema.Anything)
# slatkin_ta_ssize = Field(schema.Anything)
# entropy_ta_ssize = Field(schema.Anything)
# iqv_ta_ssize = Field(schema.Anything)
# unlabeled_config_counts_ta_ssize = Field(schema.Anything)
# num_configurations_ta_ssize = Field(schema.Anything)
# config_slatkin_ta_ssize = Field(schema.Anything)
# config_entropy_ta_ssize = Field(schema.Anything)
# config_iqv_ta_ssize = Field(schema.Anything)
# kandler_remaining_tassize = Field(schema.Anything)

def export_ta_sampled_stats():

    # METHOD NOT IMPLEMENTED YET
    return

    ## export a file with sampled statistics
    full_filename = ''
    full_filename += args.filename
    full_filename += "-tasampled-data.csv"

    sim_fields = data.mixture_model_stats.tassize_columns_to_export()

    ofile = open(full_filename, "wb")
    writer = csv.DictWriter(ofile, fieldnames=sim_fields, quotechar='"', quoting=csv.QUOTE_ALL)
    headers = dict((n, n) for n in sim_fields)
    writer.writerow(headers)
    cursor = data.MixtureModelStats.m.find(dict(), dict(timeout=False))

    for sample in cursor:
        pass

    ofile.close()


def export_slatkin_locus_values():
    # ## Export a full population census statistics file ###
    full_filename = ''
    full_filename += args.filename
    full_filename += "-pop-slatkin-locus-data.csv"
    pop_fields = []

    # adjust the fields for the new summary statistics
    pop_fields.append('simulation_run_id')
    pop_fields.append('model_class_label')
    pop_fields.append('innovation_rate')
    pop_fields.append('slatkin_locus_value')

    ofile = open(full_filename, "wb")
    writer = csv.DictWriter(ofile, fieldnames=pop_fields, quotechar='"', quoting=csv.QUOTE_ALL)
    headers = dict((n, n) for n in pop_fields)
    writer.writerow(headers)

    cursor = data.MixtureModelStats.m.find(dict(), dict(timeout=False))
    for sample in cursor:
        # slatkin exact
        slatkin_values = sample['slatkin_exact']

        for value in slatkin_values:
            row = dict()
            row['simulation_run_id'] = sample['simulation_run_id']
            row['model_class_label'] = sample['model_class_label']
            row['innovation_rate'] = sample['innovation_rate']
            row['slatkin_locus_value'] = value

            #log.info("sim data row: %s", row)
            writer.writerow(row)

    ofile.close()

def export_richness_locus_values():
    # ## Export a full population census statistics file ###
    full_filename = ''
    full_filename += args.filename
    full_filename += "-pop-richness-locus-data.csv"
    pop_fields = []

    # adjust the fields for the new summary statistics
    pop_fields.append('simulation_run_id')
    pop_fields.append('model_class_label')
    pop_fields.append('innovation_rate')
    pop_fields.append('richness_locus_value')

    ofile = open(full_filename, "wb")
    writer = csv.DictWriter(ofile, fieldnames=pop_fields, quotechar='"', quoting=csv.QUOTE_ALL)
    headers = dict((n, n) for n in pop_fields)
    writer.writerow(headers)

    cursor = data.MixtureModelStats.m.find(dict(), dict(timeout=False))
    for sample in cursor:
        # slatkin exact
        richness_values = sample['pop_richness']

        for value in richness_values:
            row = dict()
            row['simulation_run_id'] = sample['simulation_run_id']
            row['model_class_label'] = sample['model_class_label']
            row['innovation_rate'] = sample['innovation_rate']
            row['richness_locus_value'] = value

            #log.info("sim data row: %s", row)
            writer.writerow(row)

    ofile.close()





############################################################################

if __name__ == "__main__":
    setup()
    export_simulation_record()
    export_population_stats()
    export_sampled_stats()
    export_ta_sampled_stats()

    export_slatkin_locus_values()
    export_richness_locus_values()
















