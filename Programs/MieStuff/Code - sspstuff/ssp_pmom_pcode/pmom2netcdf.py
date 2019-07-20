# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 11:14:56 2015

@author: nesh
"""

import scipy.io as sio
Npmomdata = sio.loadmat('../IORs/iors_water_273 (Hybrid)/ssp_getpsd_T273_S331_Npmom.mat')
Npmom = Npmomdata['Npmomarray']
pmomdata = sio.loadmat('../IORs/iors_water_273 (Hybrid)/ssp_getpsd_T273_S331_pmom.mat')
pmom = pmomdata['pmomarray']

