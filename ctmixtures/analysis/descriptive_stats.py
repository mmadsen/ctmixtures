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


#################################################################################


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


#################################################################################


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

    def update(self, timestep):
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


#################################################################################

class TimeAveragedPopulationTraitAnalyzer(PopulationTraitAnalyzer):
    """
    Decorates a normal PopulationTraitAnalyzer class with additional tracking of
    time averaged trait counts, using the TimeAverager classes from the pytransmission
    package (https://github.com/mmadsen/pytransmission).

    None of the original methods are overridden except update(), which calls the wrapped
    PopulationTraitAnalyzer and then passes updated trait counts to all of the contained
    time averager classes.  All of the original statistics are available for the
    non-time averaged population observations, while in parallel a set of methods
    are exposed for accessing time averaged versions of the same stats.  The latter return
    the same data format as the population methods, but embedded in a dict where the key
    is the duration of time averaging in generations (not model ticks, which will vary
    with population size).  By convention, time averaged data are returned from the "ending"
    time averager, which would include the final tick/generation of the simulation, which
    makes the statistics comparable to the population statistics, which are calculated after
    the last call to update (typically at the final sample of the simulation run).
    """
    def __init__(self, model, starting_timeaverager, ending_timeaverager):
        # call superclass constructor
        super(TimeAveragedPopulationTraitAnalyzer, self).__init__(model)
        self.ta_trackers = []
        self.ta_trackers.append(starting_timeaverager)
        self.ta_trackers.append(ending_timeaverager)
        self.ending_ta = ending_timeaverager
        self.starting_ta = starting_timeaverager

    # decorate update from the superclass to pass its results to the list of time averaging objects
    def update(self, timestep):
        super(TimeAveragedPopulationTraitAnalyzer, self).update(timestep)

        # The TA trackers expect a map with loci as keys, and dicts as values, where the value
        # dicts are dicts of trait:count. The original pop/sample trait analyzers keep a list
        # of dicts, where locus is implicit in the list position.  TODO:  commensurate data structures in the population analyzer & tatracker
        locus = 0
        countmap = dict()
        for locus_map in self.counts:
            countmap[locus] = locus_map
            locus += 1


        # if the timestep is within the intervals of any of the timeaverager objects, we record
        # both trait counts for all loci/dimensions, and the intersected configurations/cultures/class counts
        for tatracker in self.ta_trackers:
            if tatracker.is_within_intervals(timestep):
                # pass the count map to the tracker
                tatracker.record_trait_count_sample(timestep,countmap,self.culture_counts)



    def get_ta_trait_frequencies(self):
        counts = self.ending_ta.get_counts_for_generation_intervals()
        # we use the measured, not configured population size in case we do population dynamics
        popsize = self.model.agentgraph.number_of_nodes()
        freqmap = dict()
        for interval, counts_by_locus in counts.items():
            locimap = dict()
            # the value of the dict for each locus should be a Counter
            for locus, counter in counts_by_locus.items():
                locusmap = dict()
                for trait, cnt in counter.items():
                    freq = float(cnt) / (float(popsize ** 2) * float(interval))
                    locusmap[trait] = freq
                locimap[locus] = locusmap
            freqmap[interval] = locimap
        #log.debug("ending TA frequency map: %s", freqmap)
        return freqmap


    def get_ta_trait_counts(self):
        return self.ending_ta.get_counts_for_generation_intervals()


    def get_ta_trait_richness(self):
        richness_map = {}
        counts = self.ending_ta.get_counts_for_generation_intervals()
        for interval, counts_by_locus in counts.items():
            richness_by_locus = dict()
            for locus, counter in counts_by_locus.items():
                nt = len([count for count in counter.values() if count > 0])
                #log.debug("richness for interval %s locus: %: %s", interval, locus, nt)
                richness_by_locus[locus] = nt
            richness_map[interval] = richness_by_locus
        #log.debug("richness_map: %s", richness_map)
        return richness_map


    def get_ta_trait_evenness_entropy(self):
        freqmap = self.get_ta_trait_frequencies()
        popsize = self.model.agentgraph.number_of_nodes()
        entropy_map = {}
        counts = self.ending_ta.get_counts_for_generation_intervals()
        for interval, counts_by_locus in counts.items():
            entropy_by_locus = dict()
            for locus, counter in counts_by_locus.items():
                freqlist = []
                for trait, cnt in counter.items():
                    freq = float(cnt) / (float(popsize ** 2) * float(interval))
                    freqlist.append(freq)
                    entropy_by_locus[locus] = diversity_shannon_entropy(freqlist)
            entropy_map[interval] = entropy_by_locus
        #log.debug("entropy_map: %s", entropy_map)
        return entropy_map


    def get_ta_trait_evenness_iqv(self):
        freqmap = self.get_ta_trait_frequencies()
        popsize = self.model.agentgraph.number_of_nodes()
        entropy_map = {}
        counts = self.ending_ta.get_counts_for_generation_intervals()
        for interval, counts_by_locus in counts.items():
            entropy_by_locus = dict()
            for locus, counter in counts_by_locus.items():
                freqlist = []
                for trait, cnt in counter.items():
                    freq = float(cnt) / (float(popsize ** 2) * float(interval))
                    freqlist.append(freq)
                    entropy_by_locus[locus] = diversity_iqv(freqlist)
            entropy_map[interval] = entropy_by_locus
        #log.debug("entropy_map: %s", entropy_map)
        return entropy_map

    def get_ta_slatkin_exact_probability(self):
        slatkin_map = {}
        counts = self.ending_ta.get_counts_for_generation_intervals()
        for interval, counts_by_locus in counts.items():
            slatkin_by_locus = dict()
            for locus, counter in counts_by_locus.items():
                counts = [count for count in counter.values()]
                slatkin_by_locus[locus] = slatkin_exact_test(counts)
            slatkin_map[interval] = slatkin_by_locus
        #log.debug("slatkin_map: %s", slatkin_map)
        return slatkin_map

    def get_ta_unlabeled_configuration_counts(self):
        return self.ending_ta.get_configuration_counts_for_generation_intervals()

    def get_ta_unlabeled_frequency_lists(self):
        counts = self.ending_ta.get_counts_for_generation_intervals()
        # we use the measured, not configured population size in case we do population dynamics
        popsize = self.model.agentgraph.number_of_nodes()
        freqmap = dict()
        for interval, counts_by_locus in counts.items():
            locifreq = []
            # the value of the dict for each locus should be a Counter
            for locus, counter in counts_by_locus.items():
                locival = []
                for trait, cnt in counter.items():
                    freq = float(cnt) / (float(popsize ** 2) * float(interval))
                    locival.append(freq)
                locifreq.append(sorted(locival, reverse=True))
            freqmap[interval] = locifreq
            #log.debug("ending TA unlabeled frequencies: %s", freqmap)
        return freqmap

    def get_ta_number_configurations(self):
        counts = self.ending_ta.get_configuration_counts_for_generation_intervals()
        num_config_map = dict()
        for interval, map in counts.items():
            num_config_map[interval] = len(map)
        return num_config_map


    def get_ta_unlableled_count_lists(self):
        ulc_map = dict()
        counts = self.ending_ta.get_counts_for_generation_intervals()
        for interval, counts_by_locus in counts.items():
            ulc_by_locus = dict()
            for locus, counter in counts_by_locus.items():
                counts = [count for count in counter.values()]
                ulc_by_locus[locus] = [count for count in counter.values()]
            ulc_map[interval] = ulc_by_locus
        #log.debug("unlabeled counts: %s", ulc_map)
        return ulc_map


    def get_ta_configuration_slatkin_test(self):
        slatkin_map = dict()
        counts = self.ending_ta.get_configuration_counts_for_generation_intervals()
        for interval, map in counts.items():
            slatkin_map[interval] = slatkin_exact_test(sorted(map.values(), reverse=True))
        return slatkin_map

    def get_ta_kandler_remaining_traits_per_locus(self):
        start_counts = self.starting_ta.get_counts_for_generation_intervals()
        #log.debug("start_counts: %s", start_counts)
        end_counts = self.ending_ta.get_counts_for_generation_intervals()
        #log.debug("end_counts: %s", start_counts)
        remaining_by_interval = dict()

        for interval, counts_by_locus in start_counts.items():
            remaining = []
            for locus, counter in counts_by_locus.items():
                start_counter = start_counts[interval][locus]
                end_counter = end_counts[interval][locus]
                remaining_set = start_counter & end_counter
                remaining.append(len(remaining_set))
            remaining_by_interval[interval] = remaining
        #log.debug("ta kandler remaining: %s", remaining_by_interval)
        return remaining_by_interval



    def get_ta_trait_frequencies_dbformat(self):
        pass

#################################################################################

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


    def update(self, timestep):
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


#################################################################################



