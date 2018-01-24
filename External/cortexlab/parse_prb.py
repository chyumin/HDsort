#!/usr/bin/python

import scipy.io as sio
import numpy as np
import sys

f = sys.argv[1]
#f = '/Volumes/hierlemann/intermediate_data/Mea1k/rolandd/cortexlab/set6/imecToWhisper.prb'

prb = {}
execfile(f, {}, prb)
matPrb = {}

## channel_groups
channel_groups_ = prb['channel_groups']
k = channel_groups_.keys()
N = len(k)
assert N is 1
channel_groups = channel_groups_[k[0]]

## channels:
matPrb['channels'] = np.array(channel_groups['channels'])

## graph:
matPrb['graph'] = np.array(channel_groups['graph'])

## geometry
geometry = np.empty((0, 3))
for key, value in channel_groups['geometry'].items():
	geometry = np.append(geometry, [[key, value[0], value[1]]], axis=0)

matPrb['geometry'] = geometry

sio.savemat(f, matPrb, appendmat=True)
