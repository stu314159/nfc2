import FluidChannel as fc
import numpy as np

# overall channel dimensions
aLx_p = 0.203 # meters
aLy_p = 0.253 # meters
aLz_p = 2.00 # meters
aN_divs = 141

# cavity parameters
cDepth = 0.240 # meters
cStart = 0.914 # meters
cEnd = 1.120 # meters

cavity = fc.ChannelCavity(cDepth,cStart,cEnd)
testChan = fc.FluidChannel(Lx_p = aLx_p, Ly_p = aLy_p, Lz_p = aLz_p,
                           N_divs = aN_divs, obst = cavity)
testChan.write_mat_file('ChanCavityTest')
testChan.write_bc_vtk()


