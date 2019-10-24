#!/usr/bin/env python
# -*- coding: utf-8 -*-

from base import _core
import numpy as _numpy


class Particle(_core.Particle):

    def GetFraction(self):
        data = _core.VectorDouble()
        self.CollectFraction(data)
        return _numpy.array(data, dtype = _numpy.float64)

    def Set(self, mass, fraction):
        fraction = _numpy.array(fraction, dtype = _numpy.float64)

        if len(fraction) != self.GetNs():
            raise Exception("Length of fraction list must equal the number of species.")
        if fraction.sum() != 1:
            raise Exception("Sum of fractions must be unity.")

        fraction_vect = _core.VectorDouble(self.GetNs())
        for s in range(self.GetNs()):
            fraction_vect[s] = fraction[s]

        self.SetData(mass, fraction_vect)
