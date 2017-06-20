# test_geom.py
'''
Simple geometry description to test nfc workflow

'''

import FluidChannel as fc

testChannel = fc.FluidChannel(Lx_p = 2.,Ly_p = 3., Lz_p = 14.,N_divs = 50,
    obst = fc.EmptyChannel(3.))
testChannel.write_bc_vtk()
