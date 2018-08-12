#[Rm, zm, s, tau, jms, theta_max, Iece] = ece_intensity(Rp, zp, theta, Rc, omega, m, F_B, F_Te, F_ne)
# M.J. Choi (mjchoi@nfri.re.kr)
# CC BY-NC-SA

# all mks units except Te
# Rp : R coordinates on beam path [m]
# zp : z coordinate on beam path [m]
# theta : angle between field and beam path [rad]
# omega : emission frequency [rad/s]
# m : harmonic number
# B : B field function in R-z coordinate [T]
# F_Te : 2d TriScatteredInterp Te function in R-z coordinates [J]
# F_ne : 2d TriScatteredInterp ne function in R-z coordinates [m^-3]
# @(x,y) : (R,z) coordinate

# it will calculate emission profile of frequency omega along s
# Rm : maximum emission position [m]
# zm : maximum emission position [m]
# Iece : ece intensity [W/m^2 Hz]
## all mks units, Temperature in [J]
