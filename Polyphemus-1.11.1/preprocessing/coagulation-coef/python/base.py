#!/usr/bin/env python
# -*- coding: utf-8 -*-

if "_mpi" in __file__:
    import coefficient_repartition_mpi_core as _core
else:
    import coefficient_repartition_core as _core


class CoefficientRepartitionBase(_core.CoefficientRepartitionBase):
    @staticmethod
    def Init(configuration_file, config_type = "default"):
        _core.CoefficientRepartitionBase_Init(configuration_file, config_type)


def Init(configuration_file, config_type = "default"):
    CoefficientRepartitionBase.Init(configuration_file, config_type)

