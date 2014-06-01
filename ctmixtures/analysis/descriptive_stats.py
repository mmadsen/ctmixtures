#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

import logging as log
from collections import defaultdict
import slatkin
import numpy as np
import random
import math as m
import pprint as pp




# I define these functions outside classes, so each class can reuse them and I can test them once and not dup code.

# TODO:  There is a naked constant here, and I feel bad about it.

def slatkin_exact_test(count_list):
    (prob, theta) = slatkin.montecarlo(100000, count_list, len(count_list))
    #log.debug("slatkin prob: %s  theta: %s", prob, theta)
    return prob


def diversity_shannon_entropy(freq_list):
    k = len(freq_list)
    sw = 0.0
    for i in range(0, k):
        sw += freq_list[i] * m.log(freq_list[i])
    if sw == 0:
        return 0.0
    else:
        return sw * -1.0


def diversity_iqv(freq_list):
    k = len(freq_list)

    if k <= 1:
        return 0.0

    isum = 1.0 - _sum_squares(freq_list)
    factor = float(k) / (float(k) - 1.0)
    iqv = factor * isum

    #logger.debug("k: %s  isum: %s  factor: %s  iqv:  %s", k, isum, factor, iqv)
    return iqv


def _sum_squares(freq_list):
    ss = 0.0
    for p in freq_list:
        ss += p ** 2.0
    return ss


## Analysis classes ##


class PopulationTraitAnalyzer(object):
    """
    Analyzer for trait counts and frequencies across the entire population.  At each
    call to calculate_trait_frequencies(), the analyzer looks at the state
    of the agent population and stores frequencies.  Subsequent calls to
    get methods will return frequencies, richness, or the Shannon entropy
    measure of evenness for the frequencies.

    To use this over time, call calculate_trait_frequencies() when you
    want a sample, and then the various get_* methods to return the
    desired metrics.

    """

    def __init__(self, model):
        self.model = model
        self.total_traits = model.agentgraph.number_of_nodes()

        # snapshots for calculating trait survival between two points or intervals
        self._snapshot_one = dict()
        self._snapshot_two = dict()

    def get_trait_frequencies(self):
        return self.freq

    def get_trait_counts(self):
        return self.counts

    def get_trait_frequencies_dbformat(self):
        # transform into the list of dicts that's more convenient to stuff into mongodb
        db_freq = []
        for locus in self.freq:
            l = []
            for key,val in locus.items():
                l.append(dict(trait=str(key),freq=val))
            db_freq.append(l)
        #log.debug("counts: %s", stored_counts)
        return db_freq

    def get_trait_richness(self):
        """
        Returns the number of traits with non-zero frequencies
        """
        richness = []
        for locus in self.freq:
            richness.append(len([freq for freq in locus.values() if freq > 0]))
        return richness

    def get_trait_evenness_entropy(self):
        entropy = []
        for locus in self.freq:
            entropy.append(diversity_shannon_entropy(locus.values()))
        return entropy

    def get_trait_evenness_iqv(self):
        iqv = []
        for locus in self.freq:
            iqv.append(diversity_iqv(locus.values()))
        return iqv

    def get_slatkin_exact_probability(self):
        slatkin = []
        for locus in self.counts:
            cnt = sorted(locus.values(), reverse=True)
            #log.info("cnt: %s", cnt)
            slatkin.append(slatkin_exact_test(cnt))
        return slatkin

    def get_unlabeled_frequency_lists(self):
        f = []

        for locus in self.freq:
            f.append(sorted(locus.values(), reverse=True))
        #log.debug("unlab freq: %s", f)
        return f

    def get_unlabeled_configuration_counts(self):
        counts = sorted(self.culture_counts.values(), reverse=True)
        #log.debug("configuration counts: %s", counts)

        return counts

    def get_number_configurations(self):
        return len(self.culture_counts)

    def get_unlableled_count_lists(self):
        f = []

        for locus in self.counts:
            f.append(sorted(locus.values(), reverse=True))
        #log.debug("unlab count: %s", f)
        return f

    def get_configuration_slatkin_test(self):
        config = self.get_unlabeled_configuration_counts()
        return slatkin_exact_test(config)

    def update(self):
        nf = self.model.simconfig.num_features
        self.counts = None
        self.freq = None
        self.culture_counts = None

        #spectra = dict()  # spectra will be locus as key, value will be dicts of popcount, numtraits
        self.counts = []  # counts will be locus as index to list, each list position is dict with key=trait, value=count
        self.freq = []  # frequencies will be locus as index to list, each list position is dict with key=trait, value=freq
        self.culture_counts = defaultdict(int)

        # set up data structures for counts and frequencies
        for i in xrange(0, nf):
            self.counts.append(defaultdict(int))
            self.freq.append(defaultdict(int))

        total = self.model.agentgraph.number_of_nodes()

        for agent_id in self.model.agentgraph.nodes():
            agent_traits = self.model.agentgraph.node[agent_id]['agent'].traits
            culture = self.model.get_traits_packed(agent_traits)
            self.culture_counts[culture] += 1
            for i in xrange(0, nf):
                self.counts[i][agent_traits[i]] += 1

        for i in xrange(0, nf):
            cnt = self.counts[i]
            for trait,count in cnt.items():
                self.freq[i][trait] = float(count) / float(total)


    def kandler_survival_start(self, timestep):
        """
        Records a snapshot of traits present in the population at timestep, recording the timestep
        as well.  The snapshots are recorded as a set of the keys from the updated trait counts,
        so we do not re-calculate anything from the population itself.
        """
        nf = self.model.simconfig.num_features
        self._snapshot_one_time = timestep
        for locus in xrange(0, nf):
            self._snapshot_one[locus] = set(self.counts[locus].keys())


    def kandler_survival_stop(self, timestep):
        """
        Records a second snapshot of traits present, in preparation for calculating the remaining traits over
        the interval
        """
        nf = self.model.simconfig.num_features
        self._snapshot_two_time = timestep
        for locus in xrange(0, nf):
            self._snapshot_two[locus] = set(self.counts[locus].keys())


    def get_kandler_remaining_traits_per_locus(self):
        """
        Returns a tuple containing the time interval in steps, and a list of remaining trait counts

        NOTE:  Kandler and Steele's (2013) analysis is conducted in the Wright-Fisher framework, and thus
        may not be exactly correct in the Moran dynamics.  At a minimum, it makes predictions in generations,
        not 1/N time steps, so this should be kept in mind when matching predicted values to observed ones.
        """
        nf = self.model.simconfig.num_features
        remaining_traits = []
        interval = self._snapshot_two_time - self._snapshot_one_time
        for locus in xrange(0, nf):
            remaining_set = self._snapshot_one[locus] & self._snapshot_two[locus]
            log.debug("locus: %s snapshot one: %s  snapshot two: %s ", locus, self._snapshot_one[locus], self._snapshot_two[locus])
            log.debug("locus: %s intersection: %s", locus, remaining_set)
            remaining_traits.append(len(remaining_set))

        return (interval, remaining_traits)



class SampledTraitAnalyzer(object):
    """
    Analyzer for simultaneous samples of the population, with samples being of potentially different sizes.

    These sample sizes are specified in the configuration file as an array and are not given on the command
    line, and thus do not affect the number of individual simulation runs performed.

    """
    def __init__(self, model):
        self.model = model
        self.sc = self.model.simconfig
        self.sample_sizes = self.sc.SAMPLE_SIZES_STUDIED
        self.total_traits = model.agentgraph.number_of_nodes()


    def update(self):
        """
        Updates various statistics given samples of the population.  Each statistic is represented by a
        dict where the key is the sample size, and the value is a list or other data structure, of the same type
        used in the whole-population analyzer.

        The goal here is to make as few passes as possible, given that we have three nested loops in the
        main counting operation, and at least two nested loops for everything else.

        :return: void
        """

        self.counts = defaultdict(int)
        self.freq = defaultdict(int)
        self.culture_counts = dict()
        nf = self.model.simconfig.num_features

        # set up data structures, which have several levels:  {ssize -> [list of loci]}  where
        # each entry in "list of loci" is a dict with:  trait:  count/freq
        for ssize in self.sample_sizes:
            self.counts[ssize] = []
            self.freq[ssize] = []
            self.culture_counts[ssize] = defaultdict(int)
            for i in xrange(0, nf):
                self.counts[ssize].append(defaultdict(int))
                self.freq[ssize].append(defaultdict(int))


        # take all of the samples and store them for use
        # and then process each one for counts and
        for ssize in self.sample_sizes:
            #log.debug("sampling ssize: %s", ssize)
            sample_ids = random.sample(self.model.agentgraph.nodes(), ssize)
            for id in sample_ids:
                # for each agent, first look at the multilocus configuration and count
                # then iterate over loci and count each separately
                agent_traits = self.model.get_agent_by_id(id).traits
                culture = self.model.get_traits_packed(agent_traits)
                self.culture_counts[ssize][culture] += 1
                for locus in xrange(0, nf):
                    trait = agent_traits[locus]
                    self.counts[ssize][locus][trait] += 1


        #log.info("counts for all sample sizes: %s", pp.pformat(self.counts))

        # now calculate frequencies
        for ssize in self.sample_sizes:
            for locus in xrange(0, nf):
                for trait, count in self.counts[ssize][locus].items():
                    self.freq[ssize][locus][trait] = float(count) / float(ssize)


        #log.info("freq for all sample sizes: %s", pp.pformat(self.freq))

    #### END update()


    def get_unlabeled_freq_by_ssize(self):
        res = dict()
        for ssize in self.sample_sizes:
            res[str(ssize)] = []
            for locus in self.freq[ssize]:
                vals = sorted(locus.values(),reverse=True)
                res[str(ssize)].append(vals)
        return res


    def get_unlabeled_counts_by_ssize(self):
        res = dict()
        for ssize in self.sample_sizes:
            res[str(ssize)] = []
            for locus in self.counts[ssize]:
                vals = sorted(locus.values(),reverse=True)
                res[str(ssize)].append(vals)
        return res

    def get_unlabeled_configuration_counts_by_ssize(self):
        res = dict()
        for ssize in self.sample_sizes:
            res[str(ssize)] = sorted(self.culture_counts[ssize].values(), reverse=True)
        return res

    def get_num_configurations_by_ssize(self):
        res = dict()
        for ssize in self.sample_sizes:
            res[str(ssize)] = len(self.culture_counts[ssize])
        return res

    def get_configuration_slatkin_by_ssize(self):
        res = dict()
        for ssize in self.sample_sizes:
            cnt = sorted(self.culture_counts[ssize].values(), reverse=True)
            res[str(ssize)] = slatkin_exact_test(cnt)
        return res

    def get_entropy_by_ssize(self):
        res = dict()
        for ssize in self.sample_sizes:
            res[str(ssize)] = []
            for locus in self.freq[ssize]:
                res[str(ssize)].append(diversity_shannon_entropy(locus.values()))
        return res


    def get_iqv_by_ssize(self):
        res = dict()
        for ssize in self.sample_sizes:
            res[str(ssize)] = []
            for locus in self.freq[ssize]:
                res[str(ssize)].append(diversity_iqv(locus.values()))
        return res


    def get_richness_by_ssize(self):
        res = dict()
        for ssize in self.sample_sizes:
            res[str(ssize)] = []
            for locus in self.counts[ssize]:
                res[str(ssize)].append(len([freq for freq in locus.values() if freq > 0]))
        return res

    def get_slatkin_by_ssize(self):
        res = dict()
        for ssize in self.sample_sizes:
            res[str(ssize)] = []
            for locus in self.counts[ssize]:
                cnt = sorted(locus.values(), reverse=True)
                res[str(ssize)].append(slatkin_exact_test(cnt))
        return res





