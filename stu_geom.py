import FluidChannel as fc

ocStu = fc.FluidChannel(Lx_p = 2., Ly_p = 2., Lz_p = 8.,
                        N_divs = 19,
                        obst = fc.EllipticalScourPit(1.0,4.,0.1))
ocStu.write_mat_file()
