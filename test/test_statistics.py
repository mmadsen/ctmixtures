#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

import logging as log
import unittest
import madsenlab.ctmixtures.utils as utils
import madsenlab.ctmixtures.traits as traits
import madsenlab.ctmixtures.population as pop
import madsenlab.ctmixtures.analysis as analysis
import os
import tempfile
import pprint as pp



class PopulationInitializationTest(unittest.TestCase):
    filename = "data/test.json"

    def test_trait_counting(self):
        log.info("test_trait_counting")

        config = utils.MixtureConfiguration(self.filename)
        config.popsize = 25
        config.num_features = 4
        config.num_traits = 10
        irule = config.INTERACTION_RULE_CLASS
        parsed = utils.parse_interaction_rule_map(irule)

        tf = traits.LocusAlleleTraitFactory(config)
        lf = pop.SquareLatticeFactory(config)
        p = pop.FixedTraitStructurePopulation(config,lf,tf)

        constructed = utils.construct_rule_objects(parsed,p)
        p.interaction_rules = constructed

        p.initialize_population()


        tfa = analysis.PopulationTraitAnalyzer(p)
        tfa.update()

        dbfreq = tfa.get_trait_frequencies_dbformat()
        log.info("db format freq: %s", pp.pformat(dbfreq))

        richness = tfa.get_trait_richness()
        log.info("richness: %s", pp.pformat(richness))

        entropy = tfa.get_trait_evenness_entropy()
        log.info("entropy: %s", pp.pformat(entropy))

        iqv = tfa.get_trait_evenness_iqv()
        log.info("iqv: %s", pp.pformat(iqv))


        slatkin = tfa.get_slatkin_exact_probability()
        log.info("slatkin: %s", pp.pformat(slatkin))


        ccounts = tfa.get_unlabeled_configuration_counts()
        log.info("ccount: %s", pp.pformat(ccounts))

        unlabeled = tfa.get_unlabeled_frequency_lists()
        log.info("unlabeled: %s", unlabeled)

        # can't test equality because we're assigning initial traits randomly
        # for locus in richness:
        #     self.assertEqual(config.num_traits, locus)

        self.assertEqual(len(richness), config.num_features)


    def test_trait_samples(self):
        log.info("test_trait_counting")

        config = utils.MixtureConfiguration(self.filename)
        config.popsize = 100
        config.num_features = 2
        config.num_traits = 30
        irule = config.INTERACTION_RULE_CLASS
        parsed = utils.parse_interaction_rule_map(irule)

        tf = traits.LocusAlleleTraitFactory(config)
        lf = pop.SquareLatticeFactory(config)
        p = pop.FixedTraitStructurePopulation(config,lf,tf)

        constructed = utils.construct_rule_objects(parsed,p)
        p.interaction_rules = constructed

        p.initialize_population()


        tfa = analysis.SampledTraitAnalyzer(p)
        tfa.update()

        res = tfa.get_unlabeled_freq_by_ssize()
        log.info("unlabeled freq by ssize: %s", res)

        res = tfa.get_unlabeled_counts_by_ssize()
        log.info("unlabeled counts by ssize: %s", res)

        res = tfa.get_slatkin_by_ssize()
        log.info("slatkin by ssize: %s", res)

        res = tfa.get_richness_by_ssize()
        log.info("richness by ssize: %s", res)

        res = tfa.get_entropy_by_ssize()
        log.info("entropy by ssize: %s", res)

        res = tfa.get_iqv_by_ssize()
        log.info("iqv by ssize: %s", res)

        ccounts = tfa.get_unlabeled_configuration_counts_by_ssize()
        log.info("configuration counts by ssize: %s", pp.pformat(ccounts))

        cslat = tfa.get_configuration_slatkin_by_ssize()
        log.info("configuration slatkin by ssize: %s", cslat)

        nconfig = tfa.get_num_configurations_by_ssize()
        log.info("config richness by ssize: %s", nconfig)

        self.assertTrue(True)




if __name__ == "__main__":
    log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    unittest.main()