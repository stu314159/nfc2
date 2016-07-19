# fc_test.py 
"""
 Test script for implementation of the FluidChannel.py object library

"""

import FluidChannel as fc


#openChannel = fc.FluidChannel() # basic empty default channel
#openChannel.write_bc_vtk()

#openChannel2 = fc.FluidChannel(N_divs = 50)
#openChannel2.write_bc_vtk()

#openChannel3 = fc.FluidChannel(Lx_p = 2.,Ly_p = 3., Lz_p = 14.,
#    obst = fc.EmptyChannel(3.))
#openChannel3.write_bc_vtk()

#openChannel4 = fc.FluidChannel(wallList = ['left','right','top'], N_divs = 50)
#openChannel4.write_bc_vtk()

#openChannel5  = fc.FluidChannel(Lx_p = 1., Ly_p = 1., Lz_p = 10.,
#                    N_divs = 15,
#                    obst = fc.SphereObstruction(r = 0.2, x_c = 0.5, y_c = 0.5, z_c = 5.))
#openChannel5.write_bc_vtk()

#PipeExpand requires Diam In, Diam Out
#openChannel6 = fc.FluidChannel(wallList = [],
#Lx_p = 2., Ly_p = 2., Lz_p = 16.,
#N_divs = 61,
#obst = fc.PipeExpand(0.8,1.8))
#openChannel6.write_bc_vtk()
#openChannel6.write_mat_file()

#PipeContract requires Diam In, Diam Out
#openChannel7 = fc.FluidChannel(wallList = [],
#Lx_p = 2., Ly_p = 2., Lz_p = 8.,
#N_divs = 61,
#obst = fc.PipeContract(1.8,0.8))
#openChannel7.write_bc_vtk()
#openChannel7.write_mat_file()

#WavyBed requires x_c, z_c, cyl_rad
#openChannel8 = fc.FluidChannel(wallList = ['bottom','top'],
#Lx_p = 1., Ly_p = 1., Lz_p = 6.,
#N_divs = 51,
#obst = fc.WavyBed(0.5,3.,0.1))
#openChannel8.write_bc_vtk()
#openChannel8.write_mat_file()

#PipeTurn requires pipe diam = 0.5 = Lo
openChannel9 = fc.FluidChannel(wallList = [],
Lx_p = 2., Ly_p = 5., Lz_p = 5.,
N_divs = 61,
obst = fc.PipeTurn(0.5)
#openChannel9.write_bc_vtk()
openChannel9.write_mat_file('Turn')

#PipeOut requires diam in and length in.  Assumes channel Lx = 4 and Ly = 4
#openChannel10 = fc.FluidChannel(wallList = [],
    #Lx_p = 4., Ly_p = 4., Lz_p = 9.,
    #N_divs = 41,
    #obst = fc.PipeOut(0.5,1.0))
#openChannel10.write_bc_vtk()
#openChannel10.write_mat_file()

#Butterfly requires diam.  Assumes Lx = 1.2, Ly = 1.2, Lz = 8 (z can vary).  Diam must be equal to 1.0
#openChannel11 = fc.FluidChannel(wallList = [],
#Lx_p = 1.2, Ly_p = 1.2, Lz_p = 8.,
    #N_divs = 11,
    #obst = fc.Butterfly(1.0))
#openChannel11.write_bc_vtk()
#openChannel11.write_mat_file()

#Tee pipes requires a diameter for pipe 1 and a diameter for pipe 2.  Assumes Lx = 2, Ly = 4, Lz = 8.  
#openChannel12 = fc.FluidChannel(wallList = [],
   # Lx_p = 2., Ly_p = 4., Lz_p = 8.,
    #N_divs = 61,
    #obst = fc.Tee(1.,0.5))
#openChannel12.write_bc_vtk()
#openChannel12.write_mat_file('Tee_Change')
