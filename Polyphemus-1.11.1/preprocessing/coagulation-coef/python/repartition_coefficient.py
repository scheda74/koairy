#!/usr/bin/env python
# -*- coding: utf-8 -*-

from base import _core
import numpy as _numpy


class RepartitionCoefficient(_core.CoefficientRepartition):

    def GetCoefficient(self, i):
        data = _core.VectorDouble()
        self.CollectCoefficient(i, data)
        return _numpy.array(data, dtype = _numpy.float64)


    def GetIndexFirst(self, i):
        data = _core.VectorInt()
        self.CollectIndexFirst(i, data)
        return _numpy.array(data, dtype = _numpy.int)


    def GetIndexSecond(self, i):
        data = _core.VectorInt()
        self.CollectIndexSecond(i, data)
        return _numpy.array(data, dtype = _numpy.int)    


    def GetAllCoefficient(self):
        return [self.GetCoefficient(i) for i in range(self.GetNsize())]


    def GetAllIndexFirst(self):
        return [self.GetIndexFirst(i) for i in range(self.GetNsize())]


    def GetAllIndexSecond(self):
        return [self.GetIndexSecond(i) for i in range(self.GetNsize())]
