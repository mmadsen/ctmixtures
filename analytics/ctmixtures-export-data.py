#!/usr/bin/env python

# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

import ming
import csv
import logging as log
import argparse
import ctmixtures.data as data
import ctmixtures.analysis as analysis
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

    if int(args.debug) == 1:
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
    pop_fields.append('slatkin_locus_mean')
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
    pop_fields.append('config_entropy')
    pop_fields.append('config_iqv')
    pop_fields.append('config_neiman_tf')
    pop_fields.append('neiman_tf_locus_min')
    pop_fields.append('neiman_tf_locus_max')
    pop_fields.append('neiman_tf_locus_mean')
    #pop_fields.append('powerlaw_locus_min')
    #pop_fields.append('powerlaw_locus_mean')
    #pop_fields.append('powerlaw_locus_max')
    #pop_fields.append('config_powerlaw_exponent')


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
        row['slatkin_locus_mean'] = np.average(slatkin_values)

        # shannon entropy
        entropy_list = sample['shannon_entropy']
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

        # Now calculate entropy and IQV from the trait configuration counts
        freqlist = []
        countlist = sample['trait_configuration_counts']
        total = sum(countlist)

        for count in countlist:
            freq = float(count) / float(total)
            freqlist.append(freq)

        row['config_entropy'] = analysis.diversity_shannon_entropy(freqlist)
        row['config_iqv'] = analysis.diversity_iqv(freqlist)

        # config_neiman_tf is easy since we have the frequency list

        row['config_neiman_tf'] = analysis.neiman_tf(freqlist)

        # Calculate Neiman's t_f statistic
        # Neiman fields are (1 / sum(freq^2)) - 1
        tflist = []
        num_loci = int(sample['num_features'])

        for locus in xrange(0, num_loci):
            locus_freq_list = sample['unlabeled_frequencies'][locus]
            tflist.append(analysis.neiman_tf(locus_freq_list))

        row['neiman_tf_locus_max'] = max(tflist)
        row['neiman_tf_locus_min'] = min(tflist)
        row['neiman_tf_locus_mean'] = np.average(tflist)








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
    #sim_fields.append('config_entropy_ssize')    TODO
    #sim_fields.append('config_iqv_ssize')   TODO

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
        #log.debug("ssizes: %s", ssizes)
        for ssize in ssizes:
            #log.debug("processing sample size: %s", ssize)
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

            #log.debug("sampled data row: %s", row)
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
#
#  ta interval -> locus -> sample size -> { data }
#


def export_ta_sampled_stats():

    ## export a file with sampled statistics
    full_filename = ''
    full_filename += args.filename
    full_filename += "-tasampled-data.csv"

    sim_fields = data.mixture_model_stats.tassize_columns_to_export()
    sim_fields.append('innovation_rate')
    sim_fields.append('ta_duration')
    sim_fields.append('sample_size')
    sim_fields.append('config_slatkin_ta_ssize')
    sim_fields.append('num_configurations_ta_ssize')
    sim_fields.append('config_iqv_ta_ssize')
    sim_fields.append('config_entropy_ta_ssize')
    sim_fields.append('richness_locus_min_tassize')
    sim_fields.append('richness_locus_max_tassize')
    sim_fields.append('richness_locus_mean_tassize')
    sim_fields.append('entropy_locus_min_tassize')
    sim_fields.append('entropy_locus_max_tassize')
    sim_fields.append('entropy_locus_mean_tassize')
    sim_fields.append('iqv_locus_min_tassize')
    sim_fields.append('iqv_locus_max_tassize')
    sim_fields.append('iqv_locus_mean_tassize')
    sim_fields.append('slatkin_locus_min_tassize')
    sim_fields.append('slatkin_locus_max_tassize')
    sim_fields.append('slatkin_locus_mean_tassize')
    sim_fields.append('kandler_locus_min_tassize')
    sim_fields.append('kandler_locus_max_tassize')
    sim_fields.append('kandler_locus_mean_tassize')
    sim_fields.append('config_neiman_tf_tassize')
    sim_fields.append('neiman_tf_locus_min_tassize')
    sim_fields.append('neiman_tf_locus_max_tassize')
    sim_fields.append('neiman_tf_locus_mean_tassize')


    ofile = open(full_filename, "wb")
    writer = csv.DictWriter(ofile, fieldnames=sim_fields, quotechar='"', quoting=csv.QUOTE_ALL)
    headers = dict((n, n) for n in sim_fields)
    writer.writerow(headers)
    cursor = data.MixtureModelStats.m.find(dict(), dict(timeout=False))

    for sample in cursor:
        log.debug("sample %s", sample['simulation_run_id'])
        # all of the other fields require iterating over sample sizes and TA intervals
        # TODO TA interval isn't explicitly recorded in the database, so must infer
        ssizes = sample['sample_size']
        ta_intervals = sample['config_slatkin_ta_ssize'].keys()
        num_loci = sample['num_features']

        # first handle those statistics which do not have per-locus measurements
        for tadur in ta_intervals:
            for ssize in ssizes:
                row = dict()
                row['simulation_run_id'] = sample['simulation_run_id']
                row['model_class_label'] = sample['model_class_label']
                row['innovation_rate'] = sample['innovation_rate']
                log.debug("Processing duration %s and ssize %s", tadur, ssize)
                row['ta_duration'] = tadur
                row['sample_size'] = ssize

                row['config_slatkin_ta_ssize'] = sample['config_slatkin_ta_ssize'][str(tadur)][str(ssize)]
                row['num_configurations_ta_ssize'] = sample['num_configurations_ta_ssize'][str(tadur)][str(ssize)]
                row['config_iqv_ta_ssize'] = sample['config_iqv_ta_ssize'][str(tadur)][str(ssize)]
                row['config_entropy_ta_ssize'] = sample['config_entropy_ta_ssize'][str(tadur)][str(ssize)]

                # now handle statistics that have per-locus measurements
                # The first set are screwed up and have the locus before the ssize, which makes processing them a pain in the butt

                richness = sample['richness_ta_ssize'][str(tadur)]
                richness_list = get_list_of_stats_for_locus_and_ssize(richness, ssize, num_loci)
                row['richness_locus_min_tassize'] = min(richness_list)
                row['richness_locus_max_tassize'] = max(richness_list)
                row['richness_locus_mean_tassize'] = np.average(richness_list)

                entropy = sample['entropy_ta_ssize'][str(tadur)]
                entropy_list = get_list_of_stats_for_locus_and_ssize(entropy, ssize, num_loci)
                row['entropy_locus_min_tassize'] = min(entropy_list)
                row['entropy_locus_max_tassize'] = max(entropy_list)
                row['entropy_locus_mean_tassize'] = np.average(entropy_list)

                iqv = sample['iqv_ta_ssize'][str(tadur)]
                iqv_list = get_list_of_stats_for_locus_and_ssize(iqv, ssize, num_loci)
                row['iqv_locus_min_tassize'] = min(iqv_list)
                row['iqv_locus_max_tassize'] = max(iqv_list)
                row['iqv_locus_mean_tassize'] = np.average(iqv_list)

                slatkin = sample['slatkin_ta_ssize'][str(tadur)]
                slatkin_list = get_list_of_stats_for_locus_and_ssize(slatkin, ssize, num_loci)
                row['slatkin_locus_min_tassize'] = min(slatkin_list)
                row['slatkin_locus_max_tassize'] = max(slatkin_list)
                row['slatkin_locus_mean_tassize'] = np.average(slatkin_list)
                
                # kandler is structured to have locus data last, which is the way I did sampled and it works well
                # don't do it this way again!!
                
                kandler_list = sample['kandler_remaining_tassize'][str(tadur)][str(ssize)]
                row['kandler_locus_min_tassize'] = min(kandler_list)
                row['kandler_locus_max_tassize'] = max(kandler_list)
                row['kandler_locus_mean_tassize'] = np.average(kandler_list)


                # config neiman tf
                config_count_map = sample['unlabeled_config_counts_ta_ssize'][str(tadur)][str(ssize)]
                config_count_list = config_count_map.values()
                total = sum(config_count_list)
                config_freq_list = [(float(x)/float(total)) for x in config_count_list]
                row['config_neiman_tf_tassize'] = analysis.neiman_tf(config_freq_list)

                # neiman locus min, max, mean TODO
                locus_freq = sample['unlabeled_freq_ta_ssize'][str(tadur)]
                locus_tf_list = get_locus_stats_for_ssize(locus_freq, analysis.neiman_tf, ssize, num_loci)
                row['neiman_tf_locus_min_tassize'] = min(locus_tf_list)
                row['neiman_tf_locus_max_tassize'] = max(locus_tf_list)
                row['neiman_tf_locus_mean_tassize'] = np.average(locus_tf_list)

                #log.debug("sampled data row: %s", row)
                writer.writerow(row)

    ofile.close()


def get_list_of_stats_for_locus_and_ssize(duration_map, ssize, num_loci):
    vals = []
    for locus in xrange(0, num_loci):
        vals.append(duration_map[str(locus)][str(ssize)])
    return vals

def get_locus_stats_for_ssize(locus_map, stat_function, ssize, num_loci):
    vals = []
    for locus in xrange(0, num_loci):
        locus_values = locus_map[str(locus)][str(ssize)]
        log.debug("locus_values: %s", locus_values)
        stat = stat_function(locus_values)
        log.debug("stat: %s", stat)
        vals.append(stat)
    return vals



############################################################################
# # misc exports


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
    #export_simulation_record()
    export_population_stats()
    #export_sampled_stats()
    export_ta_sampled_stats()

    #export_slatkin_locus_values()
    #export_richness_locus_values()
















