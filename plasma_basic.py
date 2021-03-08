"""
Basic plasma parameter and functions
Written by Tongnyeol Rhee
NFRI, Korea

Updated 
col_rel : collision frequency between relativistic electron and thermal electron 
eleskin : electron skin depth
col_ei  : electron ion collision frequency from F.F.Chen Appendix A
col_tokamak : collision frequency in Tokamak from wikipedia

gyroe   : gyroradius of electron
gyroi   : gyroraduis of ions
gyroe_re: gyroradius of relativistic electron
25 June 2018
Tongnyeol Rhee

colii_callen: Bug modification, temperature unit is eV not keV
10 January 2019
Tongnyeol Rhee

colii_Goldston: Goldstone collision frequency shown in his boot eq 11.24
11 Janunary 2019
Tongnyeol Rhee

Converted to Python3 
by Tongnyeol Rhee
05 Dec. 2019
"""
import numpy as np
me = 9.10938215e-31     # electron rest mass
mp = 1.672621777e-27    # proton mass
mup = 1.007276466812    # proton Isotope mass

mi = mp/mup             # Atomic mass unit
mud = 2.01410178        # Deterium Isotope mass
mut = 3.0160492         # Tritium Isotope mass 
md  = mud*mp            # Deterium mass in kg
mt  = mut*mp            # Tritium mass in kg

cv=2.99792458e8         # speed of light im m/s
eV2J = 1.602176565e-19  # Joule energy of 1 electron volt
evstat=4.8032e-10
mec2=0.510998910e6      # rest energy of electron in eV unit
mpc2=938.272046e6       # rest energy of proton in eV unit
mu0=4.e-7*np.pi         # permeability of free space
epsilon = 8.8541878176e-12 # permittivity of free space
mp2me=1.8362e3
#kstar parameter
circumf=11.309733552923255e0
r0=1.8e0
kbolt=1.3807e-23
d2r = np.pi/180.
r2d = 180./np.pi
eme = -eV2J/me

def ev2vi_n(eV,mu):
    """
    electron volt to ion velocity in m/s unit
    input:
        eV, mu (ion mass ratio to proton)
    output:
        speed of m/s unit
    """
    return np.sqrt(2.*eV*eV2J/(mu*mp))

def ev2vi_nrl(eV,mu):
    """
    electron volt to ion velocity in m/s unit
    by using NRL formulae
    input:
        eV, mu (ion mass ratio to proton)
    output:
        speed of m/s unit
    """
    return 9.79e3/np.sqrt(mu)*np.sqrt(2.*eV)

def ev2vi(eV,mu):
    """
    electron volt to ion velocity in m/s unit
    for relativistic ions
    input:
        eV, mu (ion mass ratio to proton)
    output:
        speed of m/s unit
    """    
    return cv*np.sqrt( eV*(eV+2.e0*mu*mpc2))/(eV+mu*mpc2)

def vi2ev(v,mu):
    """
    ion velocity to eV energy 
    input:
        v is speed in m/s unit
        mu is ion mass ratio to proton
    """
    return 0.5*mu*mp*v**2/eV2J

def ev2ve(eV):
    """
    electron volt to ion velocity in m/s unit
    for relativistic ions
    input:
        eV
    output:
        speed of m/s unit
    """    
    return cv*np.sqrt( eV*(eV+2.e0*mec2))/(eV+mec2)

def wce(B):
    """
    electron cyclotron frequency in rad/s unit
    input : B
    output : angular frequency
    """
    return eme*B
def B_res_e(f_ece, harm = 1.):
    """
    resonant magnetic field strength
    input : natural frequency
    output: magnetic field strength in Tesla
    """
    return me/eV2J*f_ece*np.pi*2./harm

def wci(B,mu):
    """
    Ion cyclrotron frequency in rad/s unit
    input : B, mu
    output : angular frequency
    """
    return eV2J*B/mp/mu

def fci(B,mu):
    """
    Ion cyclotron frequencyt in /s unit
    input : B,mu
    output : frequency
    """
    return wci(B,mu)/(2.*np.pi)

def fce(B):
    """
    electron cyclotron frequency in /s unit
    input : B
    output : frequency
    """
    return wce(B)/(2.*np.pi)

def fce_re(B,E):
    """
    relativistic electron cyclotron frequency in /s unit
    input : B
          : E energy in eV unit
    output : frequency
    """
    return wce(B)/(2.*np.pi)/ev2gamma(E);

def Greenwald_limit(Ip, minor_radi):
    """
    Greenwald density limit
    input: Ip plasma current in MA unit
    minor_radi: minor radius of plasma column
    output: density in 10^20 m-3 unit
    """
    return Ip/np.pi/minor_radi**2

def Alfvenv(density, mu, B):
    """
    Alfven speed defined as B/sqrt(mu0 rho)
    input: 
        density : plasma density in m^-3 unit
        mu      : mass ratio to proton
        B       : magnetic field in Tesla unit
    output:
        Alfven speed in m/s unit
    """
    return B/np.sqrt(mu0*mu*mp*density)

def AlfvenW(density,mu,B,q,R):
    """
    Alfven eigen mode frequency calculator
    input:
        density : plasma density in m^-3 unit
        mu      : mass ratio to proton
        B       : magnetic field strength in Tesla unit
        q       : local safety factor q
        R       : major radius in m unit
    output:
        print Alfven gap frequency
    """
    freq = Alfvenv(density,mu,B)/(2.*q*R)/2./np.pi;
    print( "TAE frequency is %10.3f kHz"%(freq/1.e3))
    print( "EAE frequency is %10.3f kHz"%(2.*freq/1.e3))
    print( "NAE frequency is %10.3f kHz"%(3.*freq/1.e3))

def thermal_i(mu,Ti):
    """
    Ion thermal speed defined by sqrt(Ti/mi)
    input:
        mu  : mass ratio to proton
        Ti  : ion temperature eV unit
    output:
        Ion thermal speed in m/s unit
    """
    return 9.79*1.e5/np.sqrt(mu/Ti)/1.e2

    
def lnlambda(Te, nev, debug=True):
    """
    Coulomb logarithm of electron-electron collision
    input:
        Te electron temperature of eV unit
        nev is electron density in m^-3 unig
    """
    if debug:
        return 14.9-0.5*np.log(nev/1.e20)+np.log(Te/1.e3)
    else:
        return 15.0
    
def cole(Te, nev):
    """
    collision frequency between electron and electron
    input:
        Te electron temperature of eV unit
        nev is electron density in m^-3 unit
    Ref. Destabilization of magentosonic-whistler waves by a relativistic runaway beam Physics of Plasmas 13, 062506 (2006)
    """
    e4e2me2=eV2J**(2.5)/epsilon**(2.0)/np.sqrt(me)/4.0/np.pi/2.0**(1.5)
    return e4e2me2*nev*lnlambda(Te,nev)/Te**(1.5)

def col_rel(Te, nev):
    """
    collision frequency between relativistic and thermal electron
    input:
        ne0 : electron density of m^-3 predisruption
        Te0 : electron temperature of predisruption in eV unit
    """
    return 1./(cv/ev2ve(Te))^(3.e0)*cole(Te,nev)

def eleskin(nev):
    """
    electron skin depth
    input:
        nev : electron density in m^-3 
    """
    return cv/wpe(nev);

def col_ei(Te, nev, Zev = 1.0):
    """
    electron ion collision frequency in F.F. Chen Appendix A
    input: 
        Te : eV unit
        nev : electron density in m^-3 unit
    """
    return 2.91e-12 * Zev*nev*lnlambda(Te,nev)/Te**1.5

def col_tokamak(Te, nev, R, eps, q):
    """
    collision frequency in Tokamak
    https://en.wikipedia.org/wiki/Collisionality
    input:
        Te: eV unit
        nev : density in m^_3
        R : major radius
        eps : inverse aspect-ratio
        q : safety factor
    """
    return col_ei(Te,nev)*ev2ve(Te)*np.sqrt(2.)*eps**(-1.5)*q*R;

def gyroi(E, B, mu, Zi,pitch):
    """
    Ion gyroradius in m unit
    Input:
        E : energy in eV
        B : Magnetic field strength in Tesla unit
        mu : Ion to mass ratio
        Zi : Charge in e unit
        pitch : pitchangle of v||^2 / v^2 
    """
    V = ev2vi(E, mu);
    Vperp = V*np.sqrt(1-pitch);
    return mu * mp * Vperp / Zi / eV2J / B;
    
def gyroe(E, B, pitch):
    """
    electron gyroradius in m unit
    Input:
        E : energy in eV
        B : Magnetic field strength in Tesla unit
        pitch : pitchangle of v||^2 / v^2 
    """
    V = ev2ve(E);
    Vperp = V*np.sqrt(1-pitch);
    return me * Vperp / eV2J / B;
    
def gyroe_re(E, B, pitch):
    """
    relativistic electron gyroradius in m unit
    Input:
        E : energy in eV
        B : Magnetic field strength in Tesla unit
        pitch : pitchangle of v||^2 / v^2 
    """
    V = ev2ve(E);
    Vperp = V*np.sqrt(1-pitch);
    return ev2gamma(E) * me * Vperp / eV2J / B;
 
def Debye(Te, nev):
    """
    Debye length calculator
    input:
        Te is electron temperature in eV unit
        nev is electron density in m^-3 unit
    """
    return 7.4339403820e3*np.sqrt(Te)/np.sqrt(nev)

def coulog(Te,nev):
    """
    Coulomb logarithm in Krammer paper
    input:
        Te is electron temperature in eV unit
        nev is electron density in m^-3 unit
    """
    return np.log(12.0*np.pi*nev*(Debye(Te,nev))**3)

def def_paramt():
    """
    Krammer collision default parameter
    """
    Zeff    = 1.0
    amu     = 2.0
    mf      = mp*amu
    return Zeff, amu,mf

def vcrit(Te):
    """
    V crit for Krammaer collision paper
    input:
        Te is electron temperature in eV unit
    """
    vcrit   = 3.0*np.sqrt(np.pi)/4.*(2.*eV2J/me)**(1.5)*(me/mi)*np.sqrt(Te**3.)
    return vcrit

def collnu(Te,nev):
    """
    Collision rate with electron for Krammer paper
    input:
        Te is electron temperature in eV unig
        nev is electron density in m^-3 unit
    """
    Zeff, amu, mf   = def_paramt();
    t1  = Zeff*np.sqrt(eV2J)*eV2J**2 * np.sqrt(me)/( np.sqrt(8.) * 
            np.sqrt(np.pi**3) * epsilon**2 * mf * 3.)
    return  t1*nev*coulog(Te,nev)/np.sqrt(Te**3.)

def collpi2(Te,nev,v):
    """
    Collision rate of pitch angle scattering
    input:
        Te is electron temperature in eV unig
        nev is electron density in m^-3 unit
        v is particle vecocity in m/2 unit
    """
    return vcrit(Te)/2./v**3*collnu(Te,nev)

def collnud(Te,nev,v):
    """
    Pitch angle scattering rate by Krammer paper
    input:
        Te is electron temperature in eV unig
        nev is electron density in m^-3 unit
        v is velocity in m/s
    """
    return vcrit(Te)/2./v**3*collnu(Te,nev)

def coli_nrl(nev,Ti,mu, zi=None):
    """
    ion collision rate
    Ti: eV
    nev: density in m^-3 unit
    mu: ion mass ratio to proton
    """
    if zi == None:
        zi=1.
    return 4.80e-8*nev*lnlambda(Ti,nev)/Ti**(1.5)/np.sqrt(mu)*zi**4

def colii_callen(mu,Ti,Te,nev,Zi=None):
    """
    ion ion collision frequency derived by J.Callen
    mu: ion mass ratio to proton
    Ti : ion Temperatrue by eV
    Te : electron temperature by eV
    nev : density in m^-3
    Zi: charge density to elctron charge
    """
    if Zi==None:
        Zi=1.
    return np.sqrt(me/mp/mu)*(Te/Ti)**(1.5)*Zi**2/np.sqrt(2)*cole(Te,nev);

def colii_Goldston(mu,Ti,Te,nev,Zi=None):
    """
    ion ion collision frequency presented in Goldstone boot eq 11.24
    mu: ion mass ratio to proton
    Ti : ion Temperatrue by eV
    Te : electron temperature by eV
    nev : density in m^-3
    Zi: charge density to elctron charge
    """
    if Zi is None:
        Zi = 1.;
    d1  = nev*Zi**4*eV2J**2.5*lnlambda(Te,nev)
    d2  = 12.*np.pi**(1.5)*epsilon**2*(mp*mu)**0.5*Ti**1.5
    return d1/d2



def ev2gamma(Energy):
    v = ev2ve(Energy);
    gamma = 1./np.sqrt(1-v**2/cv**2);
    return gamma;

def wpe(nev):
    """
    electron plasma oscillation frequency
    nev: electron density
    """
    return np.sqrt(nev*eV2J**2/(me*epsilon));

def wpi(nev,mu):
    """
    ion plasma oscillation frequency
    nev: electron density
    mu: ion mass ratio to proton
    """
    return np.sqrt(nev*eV2J**2/(mp*mu*epsilon));


def RpZ2xyz(R):
    X    = np.zeros(3,dtype='float')
    X[0] = R[0]*np.cos(R[1]);
    X[1] = R[0]*np.sin(R[1]);
    X[2] = R[2];
    return X;

def xyz2RpZ(X):
    R    = np.zeros(3,dtype='float')
    R[0] = np.sqrt( X[0]**2 + X[1]**2);
    R[1] = np.arctan2(X[1], X[0]);
    R[2] = X[2];
    return R;

def Vxyz2RpZ(X, VX):
    VR    = np.zeros(3,dtype='float')
    phi   = np.arctan2(X[1],X[0]);
    VR[0] = VX[0]*np.cos(phi)    + VX[1]*np.sin(phi);
    VR[1] =-VX[0]*np.sin(phi)    + VX[1]*np.cos(phi);
    VR[2] = VX[2];
    return VR;

def VRpZ2xyz2(X, VR): 
    VX    = np.zeros(3,dtype='float')
    phi   = np.arctan2(X[1],X[0]); 
    VX[0] = VR[0]*np.cos(phi)    - VR[1]*np.sin(phi); 
    VX[1] = VR[0]*np.sin(phi)    + VR[1]*np.cos(phi); 
    VX[2] = VR[2];
    return VX;

def gyration_position(B, X, vpe, rho,phase): 
    Xout    = np.zeros(3,dtype='float')
    rpe = np.zeros(3,dtype='float'); 
    Babs = np.sqrt(np.dot(B,B)); 
    if(B[1] !=0.): 
        BxBy = B[0]/B[1]; BxBy2 = BxBy*BxBy; 
        rho2 = rho / np.sqrt(1 + BxBy2); 
        rho3 = rho / np.sqrt(1 + 1./BxBy2); 
        Xo2 = X[0] +rho2; 
        Yo2 = X[1] -rho2*BxBy; 
        Xo3 = X[0] -rho2; 
        Yo3 = X[1] +rho2*BxBy; 
        Ro2 = Xo2*Xo2 + Yo2*Yo2; 
        Ro3 = Xo3*Xo3 + Yo3*Yo3; 
        if(Ro2 >= Ro3): 
            Xout[0] =Xo2 ; 
            Xout[1] =Yo2 ; 
            rpe[0] = rho2; 
            rpe[1] = -rho2*BxBy; 
        else: 
            Xout[0] =Xo3 ; 
            Xout[1] =Yo3 ;
            rpe[0] = -rho2;
            rpe[1] = rho2*BxBy;
        rpe[2] = 0.; 
    else:
        rpe[0] = 0.;
        rpe[1] = rho;
        rpe[2] = 0.;
        if(X[1] <0):
                rpe[2] *= -1.; 
    dth = -np.mod(phase,2.*np.pi); 
    Xout = ArbitraryRotate(rpe, dth,B); 
    Xout[0] += X[0]; 
    Xout[1] += X[1]; 
    Xout[2] += X[2];
    return Xout;

def gyration_velocity( B, X, vpe, rho,  phase):
    Xout    = np.zeros(3,dtype='float')
    rpe = np.zeros(3,dtype='float')
    rc = np.zeros(3,dtype='float')
    Babs = np.sqrt(np.dot(B,B));
    if(B[1] !=0.): 
            BxBy = B[0]/B[1];
            BxBy2 = BxBy*BxBy;
            rho2 = rho / np.sqrt(1 + BxBy2);
            rho3 = rho / np.sqrt(1 + 1./BxBy2);
            Xo2 = X[0] +rho2;
            Yo2 = X[1] -rho2*BxBy;
            Xo3 = X[0] -rho2;
            Yo3 = X[1] +rho2*BxBy;
            Ro2 = Xo2*Xo2 + Yo2*Yo2;
            Ro3 = Xo3*Xo3 + Yo3*Yo3;
            if(Ro2 >= Ro3): 
                Xout[0] =Xo2 ; 
                Xout[1] =Yo2 ;
                rpe[0] = rho2;
                rpe[1] = -rho2*BxBy;
             
            else: 
                Xout[0] =Xo3 ; 
                Xout[1] =Yo3 ;
                rpe[0] = -rho2;
                rpe[1] = rho2*BxBy;
             
            rpe[2] = 0.;
     
    else:
            rpe[0] = 0.;
            rpe[1] = rho;
            rpe[2] = 0.;
            if(X[1] <0):
                 rpe[2] *= -1.;
    
    dth = -np.mod(phase,2.*np.pi); 
    rc = ArbitraryRotate(rpe, dth,B); 
    Ncoef = vpe/rho/Babs;
    Xout[0] = (rc[1]*B[2] - rc[2]*B[1])* Ncoef;
    Xout[1] = (rc[2]*B[0] - rc[0]*B[2])* Ncoef;
    Xout[2] = (rc[0]*B[1] - rc[1]*B[0])* Ncoef;
    return Xout

def ArbitraryRotate( p,  theta,  ro):
    q = np.zeros(3,dtype = 'float');
    r = np.zeros(3,dtype = 'float');
    
    rn = np.sqrt(np.dot(ro,ro));
    r  =ro/rn;
    
    costheta = np.cos(theta);
    sintheta = np.sin(theta);
    
    q[0] += (costheta + (1. - costheta) * r[0] * r[0]) * p[0];
    q[0] += ((1. - costheta) * r[0] * r[1] - r[2] * sintheta) * p[1];
    q[0] += ((1. - costheta) * r[0] * r[2] + r[1] * sintheta) * p[2];
    
    q[1] += ((1. - costheta) * r[0] * r[1] + r[2] * sintheta) * p[0];
    q[1] += (costheta + (1. - costheta) * r[1] * r[1]) * p[1];
    q[1] += ((1. - costheta) * r[1] * r[2] - r[0] * sintheta) * p[2];
    
    q[2] += ((1. - costheta) * r[0] * r[2] - r[1] * sintheta) * p[0];
    q[2] += ((1. - costheta) * r[1] * r[2] + r[0] * sintheta) * p[1];
    q[2] += (costheta + (1. - costheta) * r[2] * r[2]) * p[2];
    return q;
 
if __name__ == '__main__' :
    print("hello")
    global nu, vc
    Te  = 3.e3    #paper parameter
    nev = 5.e19   #paper parameter
    Ef  = 1.e6    #paper parameter

    Te  = 1.e2
    nev = 0.5e19
    Ef  = 100.e3
    vc  = vcrit(Te)
    nu  = collnu(Te,nev)
    print( 'Ef = 1MeV')
    print( 'Te=%7.3fe3, nev is '%(Te/1000.), nev)
    print( 'Vcrit is ',(vc)**(1./3.))
    print( 'Debye length is ',Debye(Te,nev))
    print( 'Collision nu is ',nu)
    print( 'Collision pitch angle ', collnud(Te,nev,ev2vi(Ef,2.)))


    def dXdt(t,y):
        global nu,vc
        dvdt=-nu*(y**3+vc)/y**3*y
        return dvdt

    from scipy.integrate import ode;
    r=ode(dXdt).set_integrator('vode', method='adams',with_jacobian=False,
                    max_step=10000)
    y0=ev2vi(Ef,2.)
    dt=0.00001
    r.set_initial_value( y0  , 0.)
    nt=400
    y=np.zeros(nt,dtype='float')
    t=np.zeros(nt,dtype='float')
    for i in range(nt): 
        r.integrate(r.t+dt)
        y[i]=vi2ev(r.y,2.)
        t[i]=r.t
    from pylab import *
    plot(t,y/1000.,'r-')
    xlim([0.,1.e-3])
    ylim([80.,Ef/1.e3])
    ylabel('Energy(keV)')
    xlabel('Time(s)')
    show()

