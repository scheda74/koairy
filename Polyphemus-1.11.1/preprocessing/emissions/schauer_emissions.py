#!/usr/bin/env python

import os, sys
from string   import find,ljust,replace
from re       import compile, search

p = compile('^#')

C16 = 0.0
PSOAT = 0.0
C19 = 0.0
C25 = 0.0
VOC = 0.0

f = open(sys.argv[1], 'r')

for line in f:
    line = line.replace("\n","")
    if p.search(line) == None:
        ls = line.split(":")

        kind = ls[0]
        carbons = float(ls[1])
        species = ls[2]
        gas = 0.0
        part = 0.0

        if ls[3] != "": gas = float(ls[3])
        if ls[4] != "": part = float(ls[4])

        sys.stderr.write("species " + species)

        if part == 0:
            sys.stderr.write(" is a VOC,")
            if carbons <= 13:
                sys.stderr.write(" and has carbons <= 13.\n")
                VOC += gas
            else:
                sys.stderr.write(" and has carbons > 13.\n")
        elif part > 0:
            sys.stderr.write(" is a SVOC,")
            if kind == "HAP":
                sys.stderr.write(" and is a HAP.\n")
                C16 += gas
                PSOAT += part
            elif kind == "ALK":
                sys.stderr.write(" and is a ALK,")
                if carbons <= 21:
                    sys.stderr.write(" and has carbons <= 21.\n")
                    C19 += gas
                    PSOAT += part
                else:
                    sys.stderr.write(" and has carbons > 21.\n")
                    C25 += gas
                    PSOAT += part
            elif kind == "OTH":
                 sys.stderr.write(" and is OTH.\n")
        else:
            sys.stderr.write(line)
            sys.exit(2)
f.close()

print "VOC", VOC
print "C16", C16
print "C19", C19
print "C25", C25
print "PSOAT", PSOAT
print
print "C16/VOC", C16/VOC
print "C19/VOC", C19/VOC
print "C25/VOC", C25/VOC
print "PSOAT/VOC", PSOAT/VOC
