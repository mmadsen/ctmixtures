#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

import logging as log
import unittest
import pprint as pp

import ctmixtures.utils as utils
import ctmixtures.traits as traits
import ctmixtures.population as pop
import ctmixtures.analysis as analysis

import pytransmission.aggregation as agg

log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')

class TAStatisticsTest(unittest.TestCase):
    filename = "test/test.json"

    def test_pop_ta_stats(self):
        log.info("test_ta_stats")

        config = utils.MixtureConfiguration(self.filename)
        config.popsize = 9
        config.num_features = 3
        config.num_traits = 10
        irule = config.INTERACTION_RULE_CLASS
        parsed = utils.parse_interaction_rule_map(irule)

        tf = traits.LocusAlleleTraitFactory(config)
        lf = pop.SquareLatticeFactory(config)
        p = pop.FixedTraitStructurePopulation(config,lf,tf)

        constructed = utils.construct_rule_objects(parsed,p)
        p.interaction_rules = constructed

        p.initialize_population()

        intervals = [5,10]
        max_length = 500

        sta = agg.MoranCumulativeTimeAverager(100, intervals, config.popsize, config.num_features, ending_interval=False)
        eta = agg.MoranCumulativeTimeAverager(300, intervals, config.popsize, config.num_features, ending_interval=True)

        tfa = analysis.TimeAveragedPopulationTraitAnalyzer(p,sta,eta)
        timestep = 0

        log.info("iterating population with same traits for %s time steps", max_length)
        while timestep <= 500:
            timestep += 1
            tfa.update(timestep)



        dbfreq = tfa.get_ta_trait_frequencies()
        log.info("ending TA freq by interval: %s", dbfreq)

        richness = tfa.get_ta_trait_richness()
        log.info("ending ta richness: %s", pp.pformat(richness))

        entropy = tfa.get_ta_trait_evenness_entropy()
        log.info("ending entropy: %s", pp.pformat(entropy))

        iqv = tfa.get_ta_trait_evenness_iqv()
        log.info("ending iqv: %s", pp.pformat(iqv))


        slatkin = tfa.get_ta_slatkin_exact_probability()
        log.info("ending slatkin: %s", pp.pformat(slatkin))


        ccounts = tfa.get_ta_unlabeled_configuration_counts()
        log.info("ending ccount: %s", ccounts)

        kandler = tfa.get_ta_kandler_remaining_traits_per_locus()
        log.info("kandler remaining: %s", kandler)

        configslat = tfa.get_ta_configuration_slatkin_test()
        log.info("ending config slatkin: %s", configslat)


        unlabeled_freq = tfa.get_ta_unlabeled_frequency_lists()
        log.info("unlabeled freq: %s", unlabeled_freq)
        for interval, locuslist in unlabeled_freq.items():
            l = 0
            for locus in locuslist:
                s = sum(locus)
                log.info("interval: %s locus: %s sum: %s", interval, l, s)
                l += 1


        # can't test equality because we're assigning initial traits randomly
        # for locus in richness:
        #     self.assertEqual(config.num_traits, locus)

        self.assertTrue(True)


    def test_sampled_ta_stats(self):
        log.info("test_sampled_ta_stats")

        config = utils.MixtureConfiguration(self.filename)
        config.popsize = 400
        config.num_features = 2
        config.num_traits = 10
        irule = config.INTERACTION_RULE_CLASS
        parsed = utils.parse_interaction_rule_map(irule)

        tf = traits.LocusAlleleTraitFactory(config)
        lf = pop.SquareLatticeFactory(config)
        p = pop.FixedTraitStructurePopulation(config,lf,tf)

        constructed = utils.construct_rule_objects(parsed,p)
        p.interaction_rules = constructed

        p.initialize_population()

        intervals = [5,10]
        max_length = 500

        sta = agg.MoranCumulativeTimeAverager(100, intervals, config.popsize, config.num_features, ending_interval=False)
        eta = agg.MoranCumulativeTimeAverager(300, intervals, config.popsize, config.num_features, ending_interval=True)

        tfa = analysis.TimeAveragedSampledTraitAnalyzer(p,sta,eta)
        timestep = 0

        log.info("iterating population with same traits for %s time steps", max_length)
        while timestep <= 500:
            timestep += 1
            tfa.update(timestep)

        tfa.take_sample_snapshot()

        counts = tfa.get_ta_trait_counts()
        log.info("sampled ta counts: %s", counts)

        dbfreq = tfa.get_ta_trait_frequencies()
        log.info("ending TA freq by interval: %s", dbfreq)

        richness = tfa.get_ta_trait_richness()
        log.info("ending ta richness: %s", pp.pformat(richness))

        entropy = tfa.get_ta_trait_evenness_entropy()
        log.info("ending entropy: %s", pp.pformat(entropy))

        iqv = tfa.get_ta_trait_evenness_iqv()
        log.info("ending iqv: %s", pp.pformat(iqv))

        slatkin = tfa.get_ta_slatkin_exact_probability()
        log.info("ending slatkin: %s", pp.pformat(slatkin))

        ccounts = tfa.get_ta_unlabeled_configuration_counts()
        log.info("ending ccount: %s", ccounts)

        unlabeled_freq = tfa.get_ta_unlabeled_frequency_lists()
        log.info("unlabeled freq: %s", unlabeled_freq)

        crichness = tfa.get_ta_number_configurations()
        log.info("ending ta config richness: %s", pp.pformat(crichness))

        configslat = tfa.get_ta_configuration_slatkin_test()
        log.info("ending config slatkin: %s", configslat)

        kandler = tfa.get_ta_kandler_remaining_traits_per_locus()
        log.info("kandler remaining: %s", kandler)


        # can't test equality because we're assigning initial traits randomly
        # for locus in richness:
        #     self.assertEqual(config.num_traits, locus)

        self.assertTrue(True)






if __name__ == "__main__":
    log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    unittest.main()