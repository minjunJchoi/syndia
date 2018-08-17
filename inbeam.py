from collections import OrderedDict

def write_inbeam(nmod, xf, xpoldeg, xzb, xxb):
    with open('data/inbeam_python.dat', 'w') as f:
        f.write("&edata\n")

        if nmod == 2:
            nmod = -1

        inbeam = OrderedDict()
        inbeam['xrtol'] = 3e-07 # % required rel. error
        inbeam['xatol'] = 3e-07 # % required abs. error
        inbeam['xstep'] = 1.00 # % integration step
        inbeam['npow'] = 1 # % power absorption on(1)/off(0)
        inbeam['ncd'] = 1 # % current drive calc. on(1)/off(0)
        inbeam['ianexp'] = 2 # % analytic (1) or experimental (2) equilibrium
        inbeam['ndns'] = 2 # % analytic (1) or interpolated (2) density profiles
        inbeam['nte'] = 2 # % analyt. (1) or interp. (2) electron temp. prof.
        inbeam['ncdroutine'] = 2 # % 0->Curba, 1->Lin-Liu, 2->Lin-Liu+momentum conservation
        inbeam['nshot'] = 25485 # % shot number
        inbeam['xtbeg'] = 4.500000 # % tbegin
        inbeam['xtend'] = 4.500000 # % tende
        inbeam['nmaxh'] = 5 # % maximum harmonics
        inbeam['nrela'] = 0 # % weakly relativistic (0) full relativistic (1)
        inbeam['nabsroutine'] = 0 # %
        inbeam['noout'] = 0 # % screen out (0)
        inbeam['nmod'] = nmod # % mode selection: O-mode (1), X-mode (-1)
        inbeam['xf'] = xf # % frequency om=2*pi*xf kstar
        inbeam['xtordeg'] = 0.00000 # % geom. optics injection angle
        inbeam['xpoldeg'] = -xpoldeg # % geom. optics injection angle (- : up, + : down)
        inbeam['xxb'] = xxb # % beam launching position [cm] for tracing calc. kstar
        inbeam['xyb'] = 0.00000 # % beam launching position [cm] for tracing calc. kstar
        inbeam['xzb'] = xzb # % beam launching position [cm] for tracing calc. kstar
        inbeam['xryyb'] = 10000.000 # % initial principal curvature
        inbeam['xrzzb'] = 10000.000 # % initial principal curvature
        inbeam['xwyyb'] = 2.50000 # % initial principal width kstar
        inbeam['xwzzb'] = 2.50000 # % initial principal width kstar
        inbeam['xpw0'] = 1.00000 # % initial beam power [MW]
        inbeam['xrmaj'] = 180.00000 # % major radius
        inbeam['xrmin'] = 50.00000 # % minor radius   ########## parameters below used for the analytic calculation (not necessary) #########
        inbeam['xb0'] = 1.990000 # % central toroidal magnetic field [T']
        inbeam['xdns'] = 2.7e+13 # % electron density [cm**(-3)] - core
        inbeam['edgdns'] = 1.e+12 # % electron density [cm**(-3)] - edge
        inbeam['xe1'] = 2.000000 # % exponent in the density profile - a
        inbeam['xe2'] = 1.000000 # % exponent in the density profile - b
        inbeam['xte0'] = 2.000000 # % electron temperature [keV] - core
        inbeam['xteedg'] = 0.1000000 # % electron temperature [keV] - edge
        inbeam['xe1t'] = 2.000000 # % exponent in the temperature profile - a
        inbeam['xe2t'] = 1.000000 # % exponent in the temperature profile - b
        inbeam['xdel0'] = 0.000000 # % Shafranov shift on the axis
        inbeam['xdeled'] = 0.000000 # % Shafranov shift on the edge
        inbeam['xelo0'] = 1.000000 # % elongation on the axis
        inbeam['xeloed'] = 1.000000 # % elongation on the edge kstar
        inbeam['xq0'] = 1.000000 # % safety factor on the axis
        inbeam['xqedg'] = 2.000000 # % safety factor on the edge

        for key in inbeam:
            if key == 'xqedg':
                f.write("{} = {:g}\n".format(key, inbeam[key]))
                f.write("/\n")
            else:
                f.write("{} = {:g},\n".format(key, inbeam[key]))
