import sys,os,glob
import os.path
#import otter
import re
import pandas


import lalburst, lalsimulation, lalmetaio
from pylal.antenna import response
from logmake import *

#--------------------------------------#
# Startup information.
# ------
# This should be moved to a
# configuration file post haste
#--------------------------------------#

ifos = ['H1', 'L1']
runs = ['O1']
calibration = 'C01'

lasttime = 1136762880
#lasttime = 1137213440

inj_families = ['ga']# 'sg'] # ['sg','ga'] # The types of injections for which frames must be produced.

frames = pandas.DataFrame(columns=('run', 'ifo', 'start time', 'duration'))

i = 0
for run in runs:
    for ifo in ifos:
        print run, ifo
        root = path = "/archive/frames/{}/hoft/{}".format(run, ifo)
        for r,d,f in os.walk(path):
            for file in f:
                file = file[:-4]   # strip the file extension out
                sttime = file[14:24]
                dur = file[25:]
                # Check if the frame already exists for another ifo, and add this to the list of ifos instead
                coinc = frames[(frames['start time']==sttime) &  (frames['duration']==dur)]
                if len(coinc)>0:
                    ind=coinc.index.values[0]
                    frames.ix[ind]['ifo'] = numpy.append(numpy.array(coinc['ifo'])[0], numpy.array(ifo))
                else:
                    frames.loc[i] = [run, numpy.array(ifo), sttime, dur]
                    i+= 1

frames.to_csv('frame_list.dat')
