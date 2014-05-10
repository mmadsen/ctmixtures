#!/usr/bin/env python
# Copyright (c) 2013.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""

class MoranDynamics(object):

    def __init__(self, config, model):
        self.config = config
        self.model = model
        self._timestep = 0

    @property
    def timestep(self):
        return self._timestep

    def update(self):
        """
        Implements a discrete version of a continuous time model (Moran dynamics)
        by selecting a random agent, and calling the step() method of that agent
        and incrementing the timestep of the model.
        """
        random_agent = self.model.get_random_agent()
        rule = random_agent.rule
        rule.step(random_agent, self._timestep)
        self._timestep += 1

