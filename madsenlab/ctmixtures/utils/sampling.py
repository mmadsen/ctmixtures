#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

import madsenlab.ctmixtures.analysis as analysis
import logging as log

def sample_mixture_model(model, args, config, timestep):
    tfa = analysis.PopulationTraitFrequencyAnalyzer(model)
    tfa.calculate_trait_frequencies()

    log.info("slatkin test: %s", tfa.get_slatkin_exact_probability())