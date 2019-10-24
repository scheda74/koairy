#!/usr/bin/env python
# -*- coding: utf-8 -*-

from base import _core
import numpy as _numpy
from particle import Particle as _Particle


class GeneralSection(_core.GeneralSection):

    def GetFractionIndexBoundary(self, i):
        data = _core.VectorInt()
        self.CollectFractionIndexBoundary(i, data)
        return _numpy.array(data, dtype = _numpy.int)


    def GetFractionBoundary(self, i):
        data = _core.VectorDouble()
        self.CollectFractionBoundary(i, data)
        return _numpy.array(data, dtype = _numpy.float64)


    def GetCompositionSection(self):
        composition = []
        for i in range(self.GetNd()):
            composition.append([self.GetFractionIndexBoundary(i),
                                self.GetFractionBoundary(i)])
        return composition


    def GetRandomParticle(self):
        p = _Particle()
        self.GenerateRandomParticle(p)
        return p;


    @staticmethod
    def GetDiameterDiscretization():
        data = _core.VectorDouble()
        _core.GeneralSection_CollectDiameter(data)
        return _numpy.array(data, dtype = _numpy.float64)


    @staticmethod
    def GetMassDiscretization():
        data = _core.VectorDouble()
        _core.GeneralSection_CollectMass(data)
        return _numpy.array(data, dtype = _numpy.float64)


    @staticmethod
    def GetCompositionDiscretization():
        composition = []
        for c in range(_core.GeneralSection_GetNc()):
            composition_local = []
            for i in range(_core.GeneralSection_GetNdStatic(c)):
                fraction_index = _core.VectorInt()
                _core.GeneralSection_CollectFractionIndexBoundaryStatic(c, i, fraction_index)

                fraction = _core.VectorDouble()
                _core.GeneralSection_CollectFractionBoundaryStatic(c, i, fraction)

                composition_local.append([_numpy.array(fraction_index, dtype = _numpy.int),
                                          _numpy.array(fraction, dtype = _numpy.float64)])
            composition.append(composition_local)
        return composition
