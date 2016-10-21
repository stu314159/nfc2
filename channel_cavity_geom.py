import FluidChannel as fc
import numpy as np

# overall channel dimensions
aLx_p = 5.0 # meters
aLy_p = 2.0 # meters
aLz_p = 15 # meters
aN_divs = 81

# cavity parameters
cDepth = 1.7 # meters
cStart = 9.0 # meters
cEnd = 10.7 # meters

cavity = fc.ChannelCavity(cDepth,cStart,cEnd)
testChan = fc.FluidChannel(Lx_p = aLx_p, Ly_p = aLy_p, Lz_p = aLz_p,
                           N_divs = aN_divs, obst = cavity)
testChan.write_mat_file('ChanCavityTest')
testChan.write_bc_vtk()


