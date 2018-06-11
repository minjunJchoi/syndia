function write_inbeam(nmod, xf, xpoldeg, xzb)
fileID = fopen(fullfile('C:\cygwin64\home\mjchoi\eqdsk2topfile', 'inbeam.dat'),'w');
fprintf(fileID,'&edata\n');

if nmod == 2;
    nmod = -1;
end

inbeam = {'xrtol',3e-07; ... % required rel. error
    'xatol', 3e-07; ... % required abs. error
    'xstep', 1.00; ... % integration step
    'npow', 1; ... % power absorption on(1)/off(0)
    'ncd', 1; ... % current drive calc. on(1)/off(0)
    'ianexp', 2; ... % analytic (1) or experimental (2) equilibrium
    'ndns', 2; ... % analytic (1) or interpolated (2) density profiles
    'nte', 2; ... % analyt. (1) or interp. (2) electron temp. prof.
    'ncdroutine', 2; ... % 0->Curba, 1->Lin-Liu, 2->Lin-Liu+momentum conservation
    'nshot', 25485; ... % shot number
    'xtbeg', 4.500000; ... % tbegin
    'xtend', 4.500000; ... % tende
    'nmaxh', 5; ... % maximum harmonics
    'nrela', 0; ... % weakly relativistic (0) full relativistic (1)
    'nabsroutine', 0; ... % 
    'noout', 0; ... % screen out (0)
    'nmod', nmod; ... % mode selection: O-mode (1), X-mode (-1)      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    'xf', xf; ... % wave frequency om=2*pi*xf kstar      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    'xtordeg', 0.00000; ... % geom. optics injection angle
    'xpoldeg', -xpoldeg; ... % geom. optics injection angle (- : up, + : down)    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    'xxb', 239.00000; ... % beam launching parameter for tracing calc. kstar
    'xyb', 0.00000; ... % beam launching parameter for tracing calc. kstar
    'xzb', xzb; ... % beam launching parameter for tracing calc. kstar  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    'xryyb', 10000.000; ... % initial principal curvature
    'xrzzb', 10000.000; ... % initial principal curvature
    'xwyyb', 2.50000; ... % initial principal width kstar
    'xwzzb', 2.50000; ... % initial principal width kstar
    'xpw0', 1.00000; ... % initial beam power [MW]
    'xrmaj', 180.00000; ... % major radius
    'xrmin', 50.00000; ... % minor radius
    %%%%% parameters below used for the analytic calculation (not necessary)
    'xb0', 1.990000; ... % central toroidal magnetic field [T]
    'xdns', 2.7e+13; ... % electron density [cm**(-3)] - core
    'edgdns', 1.e+12; ... % electron density [cm**(-3)] - edge
    'xe1', 2.000000; ... % exponent in the density profile - a
    'xe2', 1.000000; ... % exponent in the density profile - b
    'xte0', 2.000000; ... % electron temperature [keV] - core
    'xteedg', 0.1000000; ... % electron temperature [keV] - edge
    'xe1t', 2.000000; ... % exponent in the temperature profile - a
    'xe2t', 1.000000; ... % exponent in the temperature profile - b
    'xdel0', 0.000000; ... % Shafranov shift on the axis
    'xdeled', 0.000000; ... % Shafranov shift on the edge
    'xelo0', 1.000000; ... % elongation on the axis
    'xeloed', 1.000000; ... % elongation on the edge kstar
    'xq0', 1.000000; ... % safety factor on the axis
    'xqedg', 2.000000; ... % safety factor on the edge
    };

for i=1:length(inbeam)
    if i == length(inbeam)
        fprintf(fileID, '%s = %g\n',inbeam{i,1}, inbeam{i,2});
        fprintf(fileID, '/\n',inbeam{i,1}, inbeam{i,2});
    else
        fprintf(fileID, '%s = %g,\n',inbeam{i,1}, inbeam{i,2});
    end
end
    
fclose(fileID);

system('scp -P 2201 /home/mjchoi/eqdsk2topfile/inbeam.dat mjchoi@172.17.250.11:~/torbeam_ifortran/run_torbeam');
