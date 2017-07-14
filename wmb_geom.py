#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 08:54:21 2017

@author: sblair
"""

import FluidChannel as fc

#overall channel dimensions
aLx_p = 1.0
aLy_p = 1.0
aLz_p = 10.0
aNdivs = 10

# wall mounted brick parameters
x_c = aLx_p/2.;
z_c = (0.4)*aLz_p;
W = aLx_p/3.;
H = aLy_p/3.;
L = aLz_p/10.;

myObst = fc.WallMountedBrick(x_c,z_c,L,W,H);

myChan = fc.FluidChannel(Lx_p=aLx_p,Ly_p=aLy_p,Lz_p=aLz_p,obst=myObst,
                         N_divs=aNdivs)                         