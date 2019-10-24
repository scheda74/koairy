#!/usr/bin/env python
# -*- coding: utf-8 -*-

if "_mpi" in __file__:
    import coefficient_repartition_mpi_core as core
    from coefficient_repartition_mpi_core import VectorInt, VectorDouble, \
        vector_int, vector_str, Particle, GeneralSection
else:
    import coefficient_repartition_core as core
    from coefficient_repartition_core import VectorInt, VectorDouble, \
        vector_int, vector_str, Particle, GeneralSection


from base import Init, CoefficientRepartitionBase
from particle import Particle
from general_section import GeneralSection
from repartition_coefficient import RepartitionCoefficient

if core.MPI2 == 1:
    try:
        from mpi4py import MPI
    except Exception:
        CoefficientRepartitionBase.MPI_Init()
else:
    print "CoefficientRepartition built without MPI support."
