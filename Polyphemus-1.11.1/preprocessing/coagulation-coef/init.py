#!/usr/bin/python
# Filename : Mcnubtest_v1.py
import os, sys
import array

import coefficient_repartition as cr

with_composition = False # True or False
output_type = "netcdf" # netcdf, bin or txt

if (with_composition):
    cr.Init("config_composition.lua", "default")
else:
    cr.Init("config_without_composition.lua", "default")

cr.GeneralSection.GetCompositionDiscretization()
cr.GeneralSection.GetDiameterDiscretization()

coef = cr.RepartitionCoefficient()

# To compute all couples at at time:
coef.ComputeAll()

# At last write coefficients in a NEtCDF file:
if output_type == "netcdf":
    name = 'coefficient.nc'
    if os.path.isfile(name):
        os.remove(name)
    coef.WriteNetCDF(name)
elif output_type == "bin":
    name = 'coefficient.bin'
    if os.path.isfile(name):
        os.remove(name)
    coef.WriteBIN(name)
elif output_type == "txt":
    name = 'coefficient.txt'
    if os.path.isfile(name):
        os.remove(name)
    coef.WriteTXT(name)
