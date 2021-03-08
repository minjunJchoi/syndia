#!/usr/bin/env python
import numpy as np;
from geqdsk import Geqdsk;
from scipy import interpolate;
from scipy.integrate import ode;
from scipy.optimize import brentq;
from mpl_toolkits.mplot3d import Axes3D;
import plasma_basic as pb;
import matplotlib.pyplot as plt;


"""
Class for equilibrium magnetic quantities from geqdsk file
Written by Tongnyeol Rhee
NFRI, Korea


======License===================
Copyright 2018 Tongnyeol Rhee

geqdsk_dk is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

geqdsk_dk is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with geqdsk_dk.  If not, see <http://www.gnu.org/licenses/>.

A copy of the LGPL license is in COPYING.LESSER. Since this is based
on (and refers to) the GPL, this is included in COPYING.

=============== Functions ==============================================
def get_B_normal(self,R,Z):
    b field
    input: R, Z
    return: b vector (R, phi, Z)


def EquiB_package_for_dk(self,R,Z):
    Equilibrium values calculator for drift equation 
    input: R, Z
    output: Dictionary 
        data['B']       : Absolute B value        
        data['Bvec']    : Bvec ( Br, Bphi, Bz) 
        data['b']       : Bvec/B       
        data['curl_b']  : Curl B
        data['Grad_B']  : Grad B
        data['Grad_B_cross_b']  : Grad B x b
        data['b_dot_curl_b']    : b cdot Curl b
        data['Bcurl_b']         : B Crul b

def b_field(self,R,Z):
    normalized magetic field at a given R, Z postion
    input: R, Z
    output: b vector (R, phi, Z)

def B_abs(self,R,Z):
    B strength value
    input: R, Z
    output: B

def B2(self,R,Z):
    B**2 value
    input: R, Z
    output: B**2

def Curl_b(self,R,Z):
    Curl b  where b = Bvector / B
    input: R, Z
    output: Curl b R, Curl b phi, Curl b Z

def Curl_B(self,R,Z):
    Curl B
    input: R, Z
    output: Curl B R, Curl B phi, Curl B Z

def Grad_B(self,R,Z):
    Grad B
    input: R, Z
    output: Grad B R, Grad B phi, Grad B Z

def Grad_Bvec(self,R,Z):
    Grad Bvec
    input: R, Z

def J(self,R,Z):
    Current calculator from magnetic field
    input: R, Z
    output: Current from magnetic field

def Curvature(self,R,Z):
    Magnetic field curvature b dot nabla b
    input: R, Z
    output: Current from magnetic field

def b_X_Grad_lnB(self,R,Z):
    vector component for Grad B drift 
    b x Grad B / B
    input: R, Z
    output: vector components

def b_X_Curv(self,R,Z):
    vector component for Curvature drift
    b x kappa where kappa is curvature
    input: R, Z
    output: vector components

def e_drift(self,R,Z,E,Lambda):
    Electron drift velocity
    R : R position
    Z : Z position
    E : Energy in eV unit
    Lambda : v||^2 / v^2 
    return
        vR, vphi,vZ

def dXdt_Lorbit(self,t,y,arg1):
    dX/dt and dV/dt calculator in R,phi Z coordinator for Lorentz equation
    input  
        t: time
        y: [R, phi, Z, VR, dphi/dt, VZ] where Vphi = R * dphi/dt
        arg1: charge / mass
    return
        [dR/dt, dphi/dt, dZ/dt, dVR/dt, d^2phi/dt^2, dVZ/dt]

def dXdt(self,t,y,arg1):
    dX/dt and dV||/dt calculator in R,phi Z coordinator for drift equation
    input  
        t: time
        y: [R, phi, Z, V||] 
        arg1: [charge, mass/charge, m] where m is mass of test particle
    return
        [dR/dt, dphi/dt, dZ/dt, dV||/dt]

def dXdt_RE(self,t,y,arg1):
    dX/dt and dV||/dt calculator in R,phi Z coordinator for relativistic drift equation
        specially runaway electron
    input  
        t: time
        y: [R, phi, Z, V||, magnetic momentum] 
        arg1: [charge, mass/charge, m, mc2] where m is mass of test particle and mc2 is 
            rest mass energy
    return
        [dR/dt, dphi/dt, dZ/dt, dV||/dt]

def psin_RZ(self,psi_n, ns = 300):
    Calculateor R, Z value having psi_n on the line from (Rmaxis, Zmaxis) 
        to (maximum R_lcfs, Zmaxis)
    input 
        psi_n: normalized psi value from 0 to 1
        ns  : number of grids of the line from Axis to boundary
    return 
        the value of [R, Z] having psi_n

def dr(self,t,y):
    Field line tracing BR, BZ
    input
        t: time
        y: [R, Z]
    return
        [dR/dphi, dZ/dphi];
        

def ang(self,A,B):
    Angle between 2D vector A and B
    return angle in radian

def polyg(self, psi_n, ndphi = 100):   
    Calculation of polygon consisting flux surface
    input 
        psi_n:normalized psi value which we want to know
        ndphi: between points toroidal angle
    return 
        flux surface R, Z


def surf_aver(self,psi_n,ndphi=100):
    surface averging at a given nomalized flux surface
    surface data should be initialized before this calculation

def surf_aver_q(self,psi_n,ndphi=100):
    surface averging at a given nomalized flux surface
    surface data should be initialized before this calculation

def init_data(self,R,Z,data_RZ):
    initialzation surface data with grids
    R: X grid value
    Z: Y grid value
    data_RZ[X,Y]: data_RZ values at (X, Y).

def init_data_o(self,surf_int_func):
    initialization surface data.
    surf_int_func(R,Z) should return surface value

def BtBphi(self,R,Z):
    local pitch
    input 
        R: radial position
        Z: vertical position
    return 
        local pitch

def abs_grad_psi(self,R,Z):
    length of poloidal flux gradient
    input:
        R: radial position
        Z: vertical position

def Closed_Lines(self,qv,the0,ndphi = 1000): 
    closed field line
    inputs
        qv : qvalue of integer
        the0: initial theta of radian
        ndphi: how fine grid
    returns
        yi[0,:] : R values
        yi[1,:] : Z values
        yi[2,:] : phi values
        yi[3,:] : theta values

def transit_path(self, E, mass, charge, sigma,Lambda, psiN, Is_Toroidal=False,RZpos = None):
    Particle orbit during 1 poloidal return
    input
        E:      is energy in eV unit
        mass:   is proton ratio
        charge: is unit e ratio, proton is 1, electron is -1, alpha is 2
        sigma:  is parallel or antiparallel
        Lambda: is v||^2/v^2
        psiN:   is 0 to 1
        Is_Toroidal: True make toroidal return
    return
        [total path, total run time, dt, drift frequencies]
        for the passing particle, 1/total time is transit frequency

def transit_path_RE(self, E, sigma, Lambda, psiN, Is_Toroidal=False):
    Particle orbit of runawya electrons during 1 poloidal return
    input
        E:      is energy in eV unit
        mass:   is proton ratio
        charge: is unit e ratio, proton is 1, electron is -1, alpha is 2
        sigma:  is parallel or antiparallel
        Lambda: is v||^2/v^2
        psiN:   is 0 to 1
    return
        [total path, total time, dt, drift frequency]
        for the passing particle, 1/total time is transit frequency

def transit_Lorbit(self, E, mass, charge, sigma,Lambda, psiN, Is_Toroidal=False,RZpos = None,
        MaxT = None):
    Particle path of Lorentz orbit
    input
        E:      is energy in eV unit
        mass:   is proton ratio
        charge: is unit e ratio, proton is 1, electron is -1, alpha is 2
        sigma:  is parallel or antiparallel
        Lambda: is v||^2/v^2
        psi:    is 0 to 1
        MaxT:   is maximum T of dt unit 
                here dt is 0.01/f where f is cyclotron frequency
    return
        [total path]


==========Code change log===================================
________________________________
1. Add grad theta vector calculation 
2. Add drift and bounce/transit frequency calculation
3. Add dXdt for drift motion
17 May 2018
Tongnyeol Rhee

________________________________
1. Add e_drift for runaway electron drift velocity calculation
2. b_X_Grad_lnB(R,Z) for calculation of gradB drift 
3. b_X_Curv(R,Z); for curvature drift calculation
25 May 2018 
Tongneyol Rhee

________________________________
1.  def GradPsi(self,R,Z): return Grad_Psip
2.  def norm_vec(self,R,Z): return Normalized Grad_Psip
3.  e_drift_rel for runway electron drift velocity calculation with new way. 
3 Sep 2018
Tongnyeol Rhee


_________________________________
Gradpsi_norm(R,Z) function is added
Tongnyeol Rhee
27 Feb. 2017
NFRI, Korea

Rlast() function is added
Rlast function return R postition of last closed flux surface at Z=Zaxis

psiN_RZ(psiN) function is added
this function return R,Z_axis position of given normalized poloidal flux surface psiN

__________________________________
q_RX(R,Z) function return q at (R,Z)
Tongnyeol Rhee
18 May 2017
NFRI, Korea


________________________________
Add routines for Flux surface average
1. def dr(t,y)
        Field line tracing BR, BZ
2. def ang(A,B)
        Angle between 2D vector A and B
        return angle in radian
3. def polyg(psi_n, ndphi = 100)
        Calculation of polygon consisting flux surface
        input 
            psi_n:normalized psi value which we want to know
            ndphi: between points toroidal angle
        return flux surface R, Z

4. def surf_aver(psi_n,ndphi=100):
        surface averging at a given nomalized flux surface
        surface data should be initialized before this calculation

5. def surf_aver_q(psi_n,ndphi=100):
        surface averging at a given nomalized flux surface
        surface data should be initialized before this calculation

6. def init_data(R,Z,data_RZ):
        initialzation surface data with grids
        R: X grid value
        Z: Y grid value
        data_RZ[X,Y]: data_RZ values at (X, Y).

7. def init_data_o(surf_int_func)
        initialization surface data.
        surf_int_func(R,Z) should return surface value

8. def BtBphi(R,Z)
        local pitch

9. def abs_grad_psi(R,Z):
        length of poloidal flux gradient

10. def Closed_Lines(qv,the0,ndphi = 1000)
        closed field line
        inputs
            qv : qvalue of integer
            the0: initial theta of radian
            ndphi: how fine grid
        returns
            yi[0,:] : R values
            yi[1,:] : Z values
            yi[2,:] : phi values
            yi[3,:] : theta values
31 May 2018
Tongnyeol Rhee

________________________________
Add routines for relativistic Lorentz Equation
def dXdt_Lorbit_RE(self,t,y,arg1):
def transit_Lorbit_RE(self, E, Lambda, sigma,psiN, Is_Toroidal=False,RZpos = None,
            MaxT = None):
25 June 2018
Tongnyeol Rhee

------------------------------------
Some modification Rlast() function for large device like ITER
02 December 2019
Tongnyeol Rhee


------------------------------------
Converted to Python 3
05 December 2019
by Tongnyeol Rhee

===========Example========================================
#Example:
## initialization of geqdks
geqdsk_file = "/home/trhee/python/g011518.001750"  #geqdsk format file
g1 = geqdsk_dk(filename =geqdsk_file )

#the position where we want to know magnetic equilibrium information
R = 2.0; Z = 0.3

#Curl B
Curl_B = g1.Curl_B(R,Z)

#Magnetic field BR, BZ, Bphi
Bvec = g1.B_field(R,Z);

#Grad B 
Grad_B = g1.Grad_B(R,Z)

#flux surface contour of given normalized poloidal flux surface
RZ = g1.poly(psin)

#3D magnetic field line path of given integer safety factor q value
RZphi = g1.Closed_Lines(qv, np.pi/2., ndphi = 50)
fig1 = plt.figure(2)
ax = fig1.add_subplot(111,projection='3d');
ax.plot(yi[:,0]*np.cos(yi[:,2]),yi[:,0]*np.sin(yi[:,2]),yi[:,1])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

#Fast ion drift orbit
E       = 90.e3 #90keV ion
mu      = 2.    #ion mass ratio to proton mass
charge  = 1.    #ion charge ratio to e
Lambda  = 0.9   #v||^2/V^2
sigma   = 1.;   # 1 is parallel and -1 is anti parallel to magnetic field
psiN    = 0.8   # 

#Runaway electron drift orbit
E_RE    = 10.e6 #10MeV electron
sigma_RE= -1


yi, t, dt, Omeg = geq.transit_period(E, mu, charge, Lambda,sigma, psi)
yi2             = geq.transit_Lorbit(E, mu, charge, Lambda,sigma, psi,MaxT=100000)
yi3, t3, dt3, Omeg3 = geq.transit_path_RE(E_RE, Lambda, sigma_RE, psi)
T       = np.arange(yi.shape[0])*dt     #Time array for each position
T[-1]   = t
print( 'bounce/transit frequency is ', 1./t)
print( 'precession frequency is ', np.abs(yi[-1,1])/t/2./np.pi)

#plotting fast ion drift orbit
plt.figure(1)
plt.plot(geq.get('rlim'),geq.get('zlim'))
plt.plot(geq.get('rbbbs'),geq.get('zbbbs'))
plt.plot(geq.get('rmaxis'),geq.get('zmaxis'),'rx')
plt.plot(yi[:,0],yi[:,2],'k-')
plt.plot(yi2[:,0],yi2[:,2],'b-')

plt.figure(2)
plt.plot(geq.get('rlim'),geq.get('zlim'))
plt.plot(geq.get('rbbbs'),geq.get('zbbbs'))
plt.plot(geq.get('rmaxis'),geq.get('zmaxis'),'rx')
plt.plot(yi3[:,0],yi3[:,2],'k-')


plt.show()

"""

class geqdsk_dk(Geqdsk):
    def __init__(self,filename = None,gR0B0=False,BDS=True):
        """
        Contructor for drift kinetic equation equilibrium
        """
        if filename == None:
            self.data = {}
            Geqdsk.__init__(self);
        else:
            self.data = {}
            Geqdsk.__init__(self,filename=filename,gR0B0=gR0B0,BDS=BDS); 
            self.init_rbx()

    def __del__(self):
        pass

    def get_lcfs_psi(self):    #return last closed flux surface
      """
      Return last closed flux surface calculated by using interpolation method
      """
      from scipy import interpolate
      rbs = self.data['rbbbs'][0]
      zbs = self.data['zbbbs'][0]
      r   = self.data['r'][0]
      z   = self.data['z'][0]
      f=interpolate.interp2d(r,z,self.get('psirz'),kind='linear')
      nbbbs = self.data['nbbbs'][0]
      psis = np.zeros((nbbbs,),np.float64)
      for i in range(nbbbs): 
        psis[i] = f(rbs[i],zbs[i])
      return np.average(psis)
    
    def get_psi(self):        
      """
      Return psi interpolation on 2D
      psi_p = f(R, Z)
      """
      if not self.psi_inter:
          from scipy import interpolate
          r   = self.data['r'][0]
          z   = self.data['z'][0]
          self.f=interpolate.interp2d(r,z,self.get('psirz'),kind='quintic')
          self.psi_inter = True
      else:
          pass;
      return self.f
    def get_psi_normal(self):
        """
        Return psi_normal interpolation on 2D
        psi_p = f_normal(R,Z)
        """
        if not self.psi_inter_normal:
          from scipy import interpolate
          r   = self.data['r'][0]
          z   = self.data['z'][0]
          psi_mag = self.get('simag')
          psiw    = self.get('psiw')
          self.f_normal=interpolate.interp2d(r,z,(self.get('psirz')-psi_mag)/psiw,kind='quintic')
          self.psi_inter_normal = True
        else:
            pass;
        return self.f_normal
    
    
    
    def get_eps(self):
      """
      Return ellipticity
      """
      minr=np.min(self.data['rbbbs'][0])
      maxr=np.max(self.data['rbbbs'][0])
      minz=np.min(self.data['zbbbs'][0])
      maxz=np.max(self.data['zbbbs'][0])
      return (maxz-minz)/(maxr-minr)
    
    def get_delta(self):
      """
      Return dshape parameter
      """
      minr=np.min(self.data['rbbbs'][0])
      maxr=np.max(self.data['rbbbs'][0])
      acenter = (minr+maxr)*0.5
      a       = (maxr-minr)/2.
      
      minz=np.min(self.data['zbbbs'][0])
      c_ind = np.where(self.data['zbbbs'][0] == minz)
      cr = self.data['rbbbs'][0][c_ind]
      c =np.abs( acenter - cr)
      
      maxz=np.max(self.data['zbbbs'][0])
      c_ind = np.where(self.data['zbbbs'][0] == maxz)
      dr = self.data['rbbbs'][0][c_ind]
      d  = np.abs(acenter - dr)
    
      return (c+d)/2./a
    
    
    def get_a(self):
      """
      Return minor radius of last closed flux surface along R direction
      """
      minr=np.min(self.data['rbbbs'][0])
      maxr=np.max(self.data['rbbbs'][0])
      return (maxr-minr)/2
    
    def get_q95(self):       
      """
      Return q value of of 95% psi value
      """
      from scipy import interpolate
      psi = self.data['psi'][0]
      q   = self.data['qpsi'][0]
      f   = interpolate.interp1d(psi,q,kind='cubic')
      psi95 = np.min(psi)+(np.max(psi)-np.min(psi))*0.95
      return f(psi95)
    
    def q_inter(self, psi,der=0):
      """
      Return q value at a given normalized psi value or its derivative.
      Its derivative is dq/dpsi not dq/dpsi_normal
      optional keyword der = 0 : q value
                             1 : 1st derivative
                             2 : 2nd derivative  
      """
      return interpolate.splev(psi,self.data['q_inter'],der=der)/self.data['delpsi'][0]**der
    
    def q_RZ(self,R,Z):
        """
        Return q value at a given R and Z position
        """
        psi_N = self.f_normal(R,Z)
        return self.q_inter(psi_N);
    
    def q_root_inter(self, qval):
      """
      Return psi value at a given q value.
      """
      return interpolate.splev(qval,self.data['q_root_inter'])
    
    
    def g_inter(self, psi,der=0):
      """
      Return g value at a given normalized psi value
      Its derivative is dg/dpsi not dg/dpsi_normal
      optional keyword der = 0 : q value
                             1 : 1st derivative
                             2 : 2nd derivative  
      """
      return interpolate.splev(psi,self.data['g_inter'],der=der)/self.data['delpsi'][0]**der
    
    def B_field(self,R,Z,cyclic=False):
      """
      Return BR, BZ value at a given R, Z position.
      B   = grad psi * grad phi - g grad phi 
      """
      if not self.psi_inter:
         f = self.get_psi();
         del f;
      else:
         pass
      psi_RZ = (self.f(R,Z)[0]-self.data['simag'][0])/self.data['psiw'][0]
      if psi_RZ <= 1.00:
         BT = -self.g_inter(psi_RZ)
      else:
         BT = -self.data['bcentr'][0]*self.data['rcentr'][0]
    
      BR = -self.f(R,Z,dy=1)[0]
      BZ =  self.f(R,Z,dx=1)[0]
    
      if cyclic:
          return np.array([BR, BT, BZ])/R
      else: 
          return np.array([BR, BZ, BT])/R
    
    
    def B_theta(self,R,Z):
        """
        Return Btheta value at a given R,Z position
        """
        return 1./R * np.sqrt( self.f(R,Z,dy=1)[0]**2 + self.f(R,Z,dx=1)[0]**2)
        
    def Bxyz(self,x,y,z):
        """
        Return bx, by, and bz at the x,y,z coordinate
        """
        r = np.sqrt(x**2+y**2);
        z = z;
        Bvec  = self.B_field(r,z);
        phi   = np.arctan2(y,x);
        bx    = Bvec[0]*np.cos(phi) - Bvec[2]*np.sin(phi);
        by    = Bvec[0]*np.sin(phi) + Bvec[2]*np.cos(phi);
        bz    = Bvec[1];
        #print 'br,bphi,cos,sin',Bvec[0],Bvec[2],np.cos(phi),np.sin(phi)
        return np.array([bx,by,bz])
    
    def local_pitch(self,R,Z):
        """
        Return local pitch defined as Bt/Bp
        """
        B = self.B_field(R,Z);
        Bp    = np.sqrt(B[0]**2+B[1]**2)
        Bt    = B[2]
        return Bt/Bp
    
    def get_ext_lim(self,rat):
        """
        Return limiter shape
        """
        R0    = self.get('rcentr')
        Z0    = 0.
        rlim  = self.get('rlim')
        zlim  = self.get('zlim')
        nlim  = len(rlim)
        rlim_new  = np.zeros(nlim,dtype='float')
        zlim_new  = np.zeros(nlim,dtype='float')
        for i in range(nlim):
            l = np.sqrt((rlim[i]-R0)**2+(zlim[i]-Z0)**2);
            phi   = np.arctan2(zlim[i]-Z0,rlim[i]-R0)
            rlim_new[i]   = l*(1.+rat)*np.cos(phi)+R0;
            zlim_new[i]   = l*(1.+rat)*np.sin(phi)+Z0;
        return rlim_new,zlim_new
    def fs(self,Rs,Zs):
        """
        Return psis from interpolation
        Input: Rs ,Zs
        Output: psis
        """
        psis=[]
        if Rs.shape != Zs.shape :
            print("Rs and Zs has different dimension")
            return 0;
    
        for i in range(Rs.shape[0]): 
            psis += [self.f(Rs[i],Zs[i])[0]];
        return np.array(psis);
    def B_fields(self,Rs,Zs):
        """
        Return B field from interpolation
        input: Rs, Zs
        output: Brs, Bzs,Bphis
        """
        Brs=[]; Bzs=[];Bphis=[];
        if Rs.shape != Zs.shape :
            print("Rs and Zs has different dimension")
            return 0;
    
        for i in range(Rs.shape[0]): 
            Br, Bz, Bphi  = self.B_field(Rs[i],Zs[i]);
            Brs += [Br]; Bzs += [Bz] ; Bphis += [Bphi]
    
        return np.array(Brs), np.array(Bzs), np.array(Bphis);
    def B_FILD(self):
        """
        return B field at FILD1 position and FILD2 position
        """
        R_f1 = 2.2234; Z_f1 = 0.428;
        R_f2 = 2.241; Z_f2 = -0.377;
        B_f1 = self.B_field(R_f1, Z_f1);
        B_f2 = self.B_field(R_f2, Z_f2);
        return [B_f1, B_f2];
    def GradPsi(self,R,Z):
        GradPsiR = self.f(R,Z,dx=1)[0];
        GradPsiZ = self.f(R,Z,dy=1)[0];
        return np.array([GradPsiR,0.,GradPsiZ]);
    def norm_vec(self,R,Z):
        GradPsi = self.GradPsi(R,Z)
        return GradPsi/np.sqrt(np.dot(GradPsi,GradPsi));

    def Gradpsi_norm(self,R,Z,idebug = False):
        """
        return n \cdot Grad psi value
        """
        Bvec  = self.B_field(R,Z);
        BR = Bvec[0];
        BZ = Bvec[1];
        GradPsi_R =  self.f(R,Z,dx=1)[0];
        GradPsi_Z =  self.f(R,Z,dy=1)[0];
        B = np.sqrt(BR**2 + BZ**2);
        bR = BR/B;
        bZ = BZ/B;
        if idebug:
            print( bR, bZ, BR, BZ )
            print( GradPsi_R, GradPsi_Z)
        GradPsi = np.sqrt(GradPsi_R**2 + GradPsi_Z**2);
        GradPsi_pa = bR * GradPsi_R + bZ*GradPsi_Z;
        if idebug: 
            return GradPsi_pa,np.sqrt(GradPsi**2 - GradPsi_pa**2);
        else: 
            return np.sqrt(GradPsi**2 - GradPsi_pa**2);
      
    def in_SOL(self,R,Z):
        if self.psi_inter_normal:
            psiv = self.f_normal(R,Z)[0];
        else:
            self.get_psi_normal();
            psiv = self.f_normal(R,Z)[0];
        if(psiv <= 1. and self.data['zbbbs_min'] <= Z <= self.data['zbbbs_max']):
            return 1;
        else:
            return 0;
    def Rlast(self):
        if not self.Rlast_define:
          self.Rlast_define = True;
          rbbbs = self.data['rbbbs'][0];
          zbbbs = self.data['zbbbs'][0];
          rs = []; zs = [];
    
          for i in range(len(zbbbs)-1):
              if(rbbbs[i] > self.data['rmaxis'][0] and 
                      self.data['zmaxis'][0]-0.5 <= zbbbs[i] <= self.data['zmaxis'][0]+0.5 ):
                  rs += [rbbbs[i]];
                  zs += [zbbbs[i]];

          rs = np.array(rs)
          zs = np.array(zs)
          self.f_zs = interpolate.interp1d(zs,rs,kind='cubic')
        return self.f_zs(self.data['zmaxis'][0]);
       
    
    def psiN_RZ(self, psiN): 
        if not self.psi_inter_normal:
            self.get_psi_normal();
        def f_psiN(R):
            return self.f_normal(R,self.data['zmaxis'][0])[0]-psiN;
        Rv = brentq(f_psiN,self.data['rmaxis'][0],self.Rlast());
        return Rv;

    def get_B_normal(self,R,Z):
        """
        b field
        input: R, Z
        output: b vector (R, phi, Z)
        """
        bR, bZ, bT  = self.B_field(R,Z)
        B = np.sqrt(bR**2+bZ**2+bT**2);
        return np.array([bR,bT,bZ])/B


    def init_operators(self):
#        x=self.get_B_abs()
        return 1


    def EquiB_package_for_dk(self,R,Z):
        """
        Equilibrium values calculator for drift equation 
        input: R, Z
        output: Dictionary 
            data['B']       : Absolute B value        
            data['Bvec']    : Bvec ( Br, Bphi, Bz) 
            data['b']       : Bvec/B       
            data['curl_b']  : Curl B
            data['Grad_B']  : Grad B
            data['Grad_B_cross_b']  : Grad B x b
            data['b_dot_curl_b']    : b cdot Curl b
            data['Bcurl_b']         : B Crul b
        """
        data={}
        psi, psi_normal, psi_R, psi_RR, psi_RZ, psi_ZR, psi_Z, psi_ZZ, g, g_psi = self.B_cal_set(R,Z)

        BR  = -psi_Z/R
        BZ  =  psi_R/R
        BT  = -g/R
        D   = np.sqrt(psi_R**2+psi_Z**2+g**2)
        B   = D/R
        RB  = R*B
        B_R = -D/R**2 + ( psi_R*psi_RR + psi_Z*psi_ZR + g*psi_R*g_psi)/(R*D)
        B_Z = (psi_R*psi_RZ + psi_Z*psi_ZZ + g*psi_Z*g_psi)/(R*D)
        bphi_Z  = (g*B_Z/B - psi_Z*g_psi)/RB
        bR_Z    = (psi_Z*B_Z/B - psi_ZZ)/RB
        bZ_R    = (psi_RR - psi_R/R - psi_R*B_R/B)/RB
        Rbphi_R = (-psi_R*g_psi + g*B_R/B)/RB
        curl_bR = -bphi_Z
        curl_bphi   = bR_Z-bZ_R
        curl_bZ     = Rbphi_R

        data['B']               = B
        data['Bvec']            = np.array([BR, BT, BZ])
        data['b']               = data['Bvec']/data['B'];
        data['curl_b']          = np.array([curl_bR,curl_bphi,curl_bZ])
        data['Grad_B']          = np.array([B_R, 0., B_Z])
        data['Grad_B_cross_b']  = np.cross(data['Grad_B'],data['b']);
        data['b_dot_curl_b']    = np.dot(data['b'],data['curl_b']);
        data['Bcurl_b']         = data['B']*data['curl_b']
        return data


    def b_field(self,R,Z):
        """
        b field
        input: R, Z
        output: b vector (R, phi, Z)
        """
        B = self.B_field(R,Z)
        return np.array([B[0],B[2],B[1]])/np.sqrt(B.dot(B))
    
    def B_abs(self,R,Z):
        """
        B value
        input: R, Z
        output: B
        """
        return np.sqrt(self.B2(R,Z))

    def B2(self,R,Z):
        """
        B**2 value
        input: R, Z
        output: B**2
        """
        B=self.B_field(R,Z)
        return B.dot(B)

    def Curl_b(self,R,Z):
        """
        Curl b
        input: R, Z
        output: Curl b R, Curl b phi, Curl b Z
        """
        psi, psi_normal, psi_R, psi_RR, psi_RZ, psi_ZR, psi_Z, psi_ZZ, g, g_psi = self.B_cal_set(R,Z)
        D  = np.sqrt(psi_R**2+psi_Z**2+g**2)
        B  = D/R
        RB = R*B
        B_R = -D/R**2 + ( psi_R*psi_RR+psi_Z*psi_ZR+g*psi_R*g_psi)/(R*D)
        B_Z = (psi_R*psi_RZ+psi_Z*psi_ZZ+g*psi_Z*g_psi)/(R*D)
        bphi_Z = (g*B_Z/B-psi_Z*g_psi)/RB
        bR_Z   = (psi_Z*B_Z/B - psi_ZZ)/RB
        bZ_R   = (psi_RR-psi_R/R-psi_R*B_R/B)/RB
        Rbphi_R = (-psi_R*g_psi+g*B_R/B)/RB
        curl_bR = -bphi_Z
        curl_bphi = bR_Z-bZ_R
        curl_bZ = Rbphi_R
        return np.array([curl_bR,curl_bphi,curl_bZ])

    def Curl_B(self,R,Z):
        """
        Curl B
        input: R, Z
        output: Curl B R, Curl B phi, Curl B Z
        """
        psi, psi_normal, psi_R, psi_RR, psi_RZ, psi_ZR, psi_Z, psi_ZZ, g, g_psi = self.B_cal_set(R,Z)
        bphi_Z = -psi_Z/R*g_psi
        bR_Z = -psi_ZZ/R
        bZ_R = (psi_RR-psi_R/R)/R
        RBphi_R = -psi_R/R*g_psi
        curl_BR = -bphi_Z
        curl_Bphi = bR_Z - bZ_R
        curl_BZ = RBphi_R
        return np.array([curl_BR, curl_Bphi, curl_BZ])


    def Grad_B(self,R,Z):
        """
        Grad B
        input: R, Z
        output: Grad B R, Grad B phi, Grad B Z
        """
        psi, psi_normal, psi_R, psi_RR, psi_RZ, psi_ZR, psi_Z, psi_ZZ, g, g_psi = self.B_cal_set(R,Z)
        D  = np.sqrt(psi_R**2+psi_Z**2+g**2)
        B  = D/R
        RB = R*B
        B_R = -D/R**2 + ( psi_R*psi_RR+psi_Z*psi_ZR+g*psi_R*g_psi)/(R*D)
        B_Z = (psi_R*psi_RZ+psi_Z*psi_ZZ+g*psi_Z*g_psi)/(R*D)
        return np.array([B_R, 0., B_Z])

    def B_cal_set(self,R,Z):
        psi = self.f(R,Z)[0]
        psi_normal = (psi-self.get('simag'))/self.get('psiw')
        psi_R = self.f(R,Z,dx=1)[0]
        psi_RR = self.f(R,Z,dx=2)[0]
        psi_RZ = self.f(R,Z,dx=1,dy=1)[0]
        psi_ZR = psi_RZ
        psi_Z = self.f(R,Z,dy=1)[0]
        psi_ZZ = self.f(R,Z,dy=2)[0]
        g = self.g_inter(psi_normal)
        g_psi = self.g_inter(psi_normal,der=1)
        return psi, psi_normal, psi_R, psi_RR, psi_RZ, psi_ZR, psi_Z, psi_ZZ, g, g_psi


    def Grad_Bvec(self,R,Z):
        """
        Grad Bvec
        input: R, Z
        """
        psi, psi_normal, psi_R, psi_RR, psi_RZ, psi_ZR, psi_Z, psi_ZZ, g, g_psi = self.B_cal_set(R,Z);
        Grad_BR     = [ (psi_ZR-psi_Z/R)/R, 0., psi_ZZ/R]
        Grad_Bphi   = [ (psi_R*g_psi-g/R)/R, 0., psi_Z*g_psi/R]
        Grad_BZ     = [ (psi_RR-psi_R/R)/R, 0., psi_RZ/R]
        return np.array([Grad_BR,Grad_Bphi,Grad_BZ])

    def J(self,R,Z):
        """
        Current calculator from magnetic field
        input: R, Z
        output: Current from magnetic field
        """
        from plasma_basic import mu0
        return self.Curl_B(R,Z)/mu0
    
    def Curvature(self,R,Z):
        """
        Magnetic field curvature b dot nabla b
        input: R, Z
        output: Current from magnetic field
        """
        p, p_n, p_R, p_RR, p_RZ, p_ZR, p_Z, p_ZZ, g, g_p = self.B_cal_set(R,Z);
        Bvec = np.array([-p_Z/R, -g/R, p_R/R]);
        D  = np.sqrt(p_R**2+p_Z**2+g**2)
        B  = D/R
        RB = R*B

        B_R = -D/R**2 + ( p_R*p_RR+p_Z*p_ZR+g*p_R*g_p)/(R*D)
        B_Z = (p_R*p_RZ+p_Z*p_ZZ+g*p_Z*g_p)/(R*D)

        bR = Bvec[0]/B;
        bf = Bvec[1]/B;
        bZ = Bvec[2]/B;

        BR_R = -(R*p_RZ - p_Z)/R**2
        Bf_R = -(p_R*R*g_p-g)/R**2
        BZ_R =  (R*p_RR - p_R)/R**2
        BR_Z = -p_ZZ/R;
        Bf_Z = -p_Z*g_p/R;
        BZ_Z =  p_RZ/R;

        bR_R = (BR_R - bR*B_R)/B;
        bf_R = (Bf_R - bf*B_R)/B;
        bZ_R = (BZ_R - bZ*B_R)/B;
        bR_Z = (BR_Z - bR*B_Z)/B;
        bf_Z = (Bf_Z - bf*B_Z)/B;
        bZ_Z = (BZ_Z - bZ*B_Z)/B;

        KR = bR*bR_R + bf**2/R+bZ*bR_Z;
        Kf = bR*bf_R + bf*bR+bZ*bf_Z;
        KZ = bR*bZ_R + bZ*bZ_Z;
        return np.array([KR, Kf, KZ]);

    def b_X_Grad_lnB(self,R,Z):
        """
        vector component for Grad B drift 
        b x Grad B / B
        input: R, Z
        output: vector components
        """
        Grad_B = self.Grad_B(R,Z);
        Bvec = self.B_field(R,Z,cyclic=True);
        B    = np.sqrt(np.dot(Bvec,Bvec));
        bvec = Bvec/B;

        return np.cross(bvec,Grad_B)/B;

    def b_X_Curv(self,R,Z):
        """
        vector component for Curvature drift
        b x kappa where kappa is curvature
        input: R, Z
        output: vector components
        """
        Bvec = self.B_field(R,Z,cyclic=True);
        B    = np.sqrt(np.dot(Bvec,Bvec));
        bvec = Bvec/B;
        Curv = self.Curvature(R,Z);

        return np.cross(bvec,Curv);

    def e_drift(self,R,Z,E,Lambda):
        """
        Electron drift velocity
        R : R position
        Z : Z position
        E : Energy in eV unit
        Lambda : v||^2 / v^2 

        return 
            np.array(vR, vphi,vZ)
        """
        v    = pb.ev2ve(E);
        vpa2 = v**2 * Lambda;
        vpe2 = v**2 - vpa2;
        Bvec = self.B_field(R,Z);
        B    = np.sqrt(np.dot(Bvec,Bvec));
        fce_re = pb.fce_re(B,E);
        vdB  =0.5*vpe2*self.b_X_Grad_lnB(R,Z);
        vcv  =    vpa2*self.b_X_Curv(R,Z);
        return (vdB+vcv)/(fce_re*np.pi*2.)


        return 1;

    def e_drift_rel(self,R,Z,E,Lambda,sigma_RE):
        v    = pb.ev2ve(E);
        vpa2 = v**2 * Lambda;
        vpe2 = v**2 - vpa2;
        Bvec = self.B_field(R,Z);
        B    = np.sqrt(np.dot(Bvec,Bvec));
        gamma = pb.ev2gamma(E)

        vpa = sigma_RE*np.sqrt(vpa2)*gamma
        mum = pb.me*gamma**2*vpe2/2./B
        Zhe = -pb.eV2J;
        mq  = 1./pb.eme
        data = self.EquiB_package_for_dk(R,Z)

        mu_ga   = mum/gamma       #mu/gamma
        rhopa = vpa*mq/data['B']  #V_||/Omega
        star = (1.+rhopa*data['b_dot_curl_b'])
        Bstar = data['B']*star
        #dX/dT
        vpa_star = vpa / gamma / star * ( rhopa*data['curl_b'])
        vB = -mu_ga/Zhe/Bstar*data['Grad_B_cross_b']
        vpa = vpa/gamma/star*data['b']
        dxdt = vpa_star + vB

        return np.array([dxdt[0], dxdt[1]/R,dxdt[2]]), vpa


    def sample_density(self,rho,a=2.,b=2.2):
        """
        sample density profile 
        rho is normalized poloidal flux
        """
        if 1.>=rho >=0.:
            return (-a*rho**2+b)*1.e13
        else:
            return 1.e10

    def equi_p(self,x,y,z):
        R = np.sqrt(x**2+y**2);
        Z = z;
        rmin = self.get('r').min()
        rmax = self.get('r').max()
        zmin = self.get('z').min()
        zmax = self.get('z').max()

        if rmin<R<rmax and zmin<Z<zmax:
            psi = (self.f(R,Z)-self.get('simag'))/self.get('psiw');
            Ne = self.sample_density(psi);
            Te = 0.
            Zeff = 0.
            return Ne,Te,Zeff
        else:
            return 1.e10,0.,0.

    def Grad_Theta(self,R,Z,dt=None):
        if(dt == None): 
            dt = 1.e-4;
        Rax = self.data['rmaxis'][0];
        Zax = self.data['zmaxis'][0];
        dtR1 = np.arctan2(Z-Zax,R-Rax)
        dtR2 = np.arctan2(Z-Zax,R+dt-Rax)
        dtZ1 = np.arctan2(Z-Zax,R-Rax)
        dtZ2 = np.arctan2(Z+dt-Zax,R-Rax)
        dtR = dtR2 - dtR1
        dtZ = dtZ2 - dtZ1
        return np.array([dtR/dt,0.,dtZ/dt])


    def dXdt_Lorbit(self,t,y,arg1):
        qm      = arg1[0]
        Bvec    = self.B_field(y[0],y[2],cyclic=True);

        R       = y[0];
        vR      = y[3];
        vphi    = R*y[4];
        vZ      = y[5];
        vBR     = vphi*Bvec[2] - vZ  *Bvec[1];
        vBphi   = vZ  *Bvec[0] - vR  *Bvec[2];
        vBZ     = vR  *Bvec[1] - vphi*Bvec[0];

        dpdt    = y[4];

        f0  = vR;
        f1  = dpdt;
        f2  = vZ;
        f3  = R*dpdt**2 + qm*vBR;
        f4  = (-2.*vR*dpdt+qm*vBphi)/R;
        f5  = qm*vBZ;

        return np.array([f0,f1,f2,f3,f4,f5]);

    def dXdt_Lorbit_RE_cylinder(self,t,y,arg1):
        qm      = arg1[0]
        gamma   = arg1[1];
        Bvec    = self.B_field(y[0],y[2],cyclic=True);

        R       = y[0];
        vR      = y[3]/gamma;
        vphi    = R*y[4]/gamma;
        vZ      = y[5]/gamma;
        vBR     = vphi*Bvec[2] - vZ  *Bvec[1];
        vBphi   = vZ  *Bvec[0] - vR  *Bvec[2];
        vBZ     = vR  *Bvec[1] - vphi*Bvec[0];

        dpdt    = y[4]/gamma;

        f0  = vR;
        f1  = dpdt;
        f2  = vZ;
        f3  = R*gamma*dpdt**2 + qm*vBR;
        f4  = (-2.*y[3]*dpdt+qm*vBphi)/R;
        f5  = qm*vBZ;

        return np.array([f0,f1,f2,f3,f4,f5]);

    def dXdt_Lorbit_RE(self,t,y,arg1):
        qm      = arg1[0]
        gamma   = arg1[1];
        Bvec    = self.Bxyz(y[0],y[1],y[2]);

        vx      = y[3]/gamma;
        vy      = y[4]/gamma;
        vz      = y[5]/gamma;
        vBx     = vy  *Bvec[2] - vz  *Bvec[1];
        vBy     = vz  *Bvec[0] - vx  *Bvec[2];
        vBz     = vx  *Bvec[1] - vy  *Bvec[0];

        f0  = vx;
        f1  = vy;
        f2  = vz;
        f3  = qm*vBx
        f4  = qm*vBy
        f5  = qm*vBz

        return np.array([f0,f1,f2,f3,f4,f5]);

    def dXdt(self,t,y,arg1):
        vpa = y[3];
        mum = y[4];
        Zhe = arg1[0];
        mq  = arg1[1];
        mh  = arg1[2];
        data = self.EquiB_package_for_dk(y[0],y[2])
        rhopa = vpa*mq/data['B']
        star = (1.+rhopa*data['b_dot_curl_b'])
        Bstar = data['B']*star
        #dX/dT
        vpa_star = vpa / star * ( data['b']+rhopa*data['curl_b'])
        vB = -mum/Zhe/Bstar*data['Grad_B_cross_b']
        dxdt = vpa_star + vB
        #parallel velocity
        dvdt = -1./star/mh*np.dot(data['b']+rhopa*data['curl_b'],mum*data['Grad_B'])
        return np.array([dxdt[0], dxdt[1]/y[0],dxdt[2],dvdt,0.])

    def dXdt_RE(self,t,y,arg1):
        vpa = y[3];
        mum = y[4];
        Zhe = arg1[0];
        mq  = arg1[1];
        mh  = arg1[2];
        mc2 = arg1[3];
        data = self.EquiB_package_for_dk(y[0],y[2])
        gamma   = 1./np.sqrt(1.+(2.*mum*data['B']+mh*vpa**2.)/mc2)
        mu_ga   = mum*gamma       #mu/gamma
        rhopa = vpa*mq/data['B']  #V_||/Omega
        star = (1.+rhopa*data['b_dot_curl_b'])
        Bstar = data['B']*star
        #dX/dT
        vpa_star = vpa * gamma / star * ( data['b']+rhopa*data['curl_b'])
        vB = -mu_ga/Zhe/Bstar*data['Grad_B_cross_b']
        dxdt = vpa_star + vB
        #parallel velocity
        dvdt = -1./star/mh*np.dot(data['b']+rhopa*data['curl_b'],mu_ga*data['Grad_B'])

        return np.array([dxdt[0], dxdt[1]/y[0],dxdt[2],dvdt,0.])

    def psin_RZ(self,psi_n, sigma = 1.):
        """
        Calculateor R, Z value having psi_n on the line from (Rmaxis, Zmaxis) 
            to (maximum R_lcfs, Zmaxis)
        input 
            psi_n: normalized psi value from 0 to 1
            ns  : number of grids of the line from Axis to boundary
        return 
            the value of [R, Z] having psi_n
        """
        if psi_n <=0.:
            print("Please input 0<psi_n<1")
            return 0.;
        """ OLD version
        rbbbs = self.get('rbbbs')   #lcfs R
        zbbbs = self.get('zbbbs')   #lcfs Z
        f = self.get_psi_normal();  #normalized poloidal flux
        R0 = self.get('rmaxis');    #magnetic axis R
        Z0 = self.get('zmaxis');    #magnetic axis Z
        Rmax = self.get('rbbbs').max();  #maximum R

        r = np.linspace(R0,Rmax, ns); 
        r_psi = np.zeros(ns);
        for i in range(ns):
            r_psi[i] = f(r[i],Z0);
        r_inter_psi = interpolate.interp1d(r, r_psi - psi_n, kind = 'cubic')
        try: 
            rv = brentq(r_inter_psi, R0, Rmax);
        except:
            print( "Error finding flux surface")
            print( "If you give close to 0 or 1 then please use diff. psi values")
            return 0.;
        else:
            return [rv,Z0]
        """
        #### New method """
        if sigma >= 0:
               return self.f_RZ_Nr(psi_n).item(), self.get('zmaxis');
        else:
               return self.f_RZ_Nl(psi_n).item(),self.get('zmaxis');


    def dr(self,t,y):
        """
        Field line tracing BR, BZ
        """
        BR, BZ, BT = self.B_field(y[0],y[1]);
        return np.array([BR, BZ])*y[0]/BT;

    def ang(self,A,B):
        """
        Angle between 2D vector A and B
        return angle in radian
        """
        a = np.sqrt(A[0]**2 + A[1]**2);
        b = np.sqrt(B[0]**2 + B[1]**2);
        ab = A[0]*B[0]+A[1]*B[1];
        return np.arccos(ab/(a*b));


    def polyg(self, psi_n, ndphi = 100):   
        """
        Calculation of polygon consisting flux surface
        input 
            psi_n:normalized psi value which we want to know
            ndphi: between points toroidal angle
        return flux surface R, Z

        """
        if psi_n <=0.:
            print("Input 0 < psi_n < 1 please")
            return 0.;
        Rax, Zax = self.get('rmaxis'), self.get('zmaxis');
        Ro = self.psiN_RZ(psi_n)
        Zo = Zax;
        dphi = np.pi*2./float(ndphi)
        theta_rot = 0.;

        r = ode(self.dr).set_integrator('vode',method = 'adams',with_jacobian = False, 
                max_step = 100000);

        r.set_initial_value([Ro,Zo],0.);

        R0, Z0 = Ro, Zo;
        yi = np.array([R0,Z0,0.]);
        while theta_rot < 2.*np.pi:
            r.integrate(r.t+dphi);
            R1, Z1 = r.y[0], r.y[1];
            dtheta = self.ang([R0-Rax,Z0-Zax],[R1-Rax,Z1-Zax]);
            theta_rot += np.abs(dtheta);
            R0, Z0 = R1, Z1;
            yi = np.vstack((yi,np.array([r.y[0],r.y[1],r.t])));
        if theta_rot == 2.*np.pi: 
            return yi;
        else:
            return np.vstack((yi[0:-1,:],np.array([Ro,Zo,r.t])));
    
    def polyg_q(self,qv,ndphi=100):
        """
        Calcuation of polygon consisting flux surface of given q value
        input
            qv      : q value 
            ndphi   : between points toroidal angle
        return
            flux surface R,Z
        """
        if self.data['qpsi'][0][0] < qv < self.data['qpsi'][0][-1]: 
            psin = self.q_root_inter(qv);
            return self.polyg(psin,ndphi)
        else:
            print(" q value is not in the equilibrium range")
            return 0;

    def surf_aver(self,psi_n,ndphi=100):
        """
        surface averging at a given nomalized flux surface
        surface data should be initialized before this calculation
        """
        ys = self.polyg(psi_n,ndphi=ndphi);
        ns = ys.shape[0];
        yc = np.zeros([ns-1,2],dtype =float)
        lc = np.zeros(ns-1,dtype = float)
        ns1 = ns - 1;
        aver = 0.;
        for i in range(ns1):
            yc[i,:] = (ys[i,0] + ys[i+1,0])*0.5, (ys[i,1]+ys[i+1,1])*0.5 #center R,Z
            lc[i] = np.sqrt((ys[i,0] - ys[i+1,0])**2 + (ys[i,1]-ys[i+1,1])**2) #length
            aver += self.surf_data(yc[i,0],yc[i,1])*lc[i]

        return aver/lc.sum();
    def surf_aver_q(self,psi_n,ndphi=100):
        """
        surface averging at a given nomalized flux surface
        surface data should be initialized before this calculation
        """
        ys = self.polyg(psi_n,ndphi=ndphi);
        ns = ys.shape[0];
        yc = np.zeros([ns-1,2],dtype = float)
        lc = np.zeros(ns-1,dtype = float)
        ns1 = ns - 1;
        aver = 0.;
        length = 0.;
        for i in range(ns1):
            yc[i,:] = (ys[i,0] + ys[i+1,0])*0.5, (ys[i,1]+ys[i+1,1])*0.5 #center R,Z
            dl = np.sqrt((ys[i,0]-ys[i+1,0])**2 + (ys[i,1]-ys[i+1,1])**2)
            length += dl;
            aver += np.abs(self.surf_data(yc[i,0],yc[i,1]))*dl

        return aver /np.pi/2.

    def init_data(self,R,Z,data_RZ):
        """
        initialzation surface data with grids
        R: X grid value
        Z: Y grid value
        data_RZ[X,Y]: data_RZ values at (X, Y).
        """
        self.surf_data=interpolate.interp2d(R,Z,RZ.transpose,kind='linear')

    def init_data_o(self,surf_int_func):
        """
        initialization surface data.
        surf_int_func(R,Z) should return surface value
        """
        self.surf_data = surf_int_func;
    def BtBphi(self,R,Z):
        """
        local pitch
        """
        BR, BZ, BT = self.B_field(R,Z);
        return BT/R/np.sqrt(BR**2+BZ**2);


    def abs_grad_psi(self,R,Z):
        """
        length of poloidal flux gradient
        """
        dpsi_R = self.psi(R,Z,dx=1)[0];
        dpsi_Z = self.psi(R,Z,dy=1)[0];
        return np.sqrt(dpsi_R**2 + dpsi_Z**2);

    def Closed_Lines(self,qv,the0,ndphi = 1000): 
        """
        closed field line
        inputs
            qv : qvalue of integer
            the0: initial theta of radian
            ndphi: how fine grid
        returns
            yi[0,:] : R values
            yi[1,:] : Z values
            yi[2,:] : phi values
            yi[3,:] : theta values
         """
        if( qv > self.get('qpsi').max() or  qv < self.get('qpsi').min()):
            print("please check q value minimum and maximum")
            return 1
        else:
            pass
        psiv = self.q_root_inter(qv) 
        R0 = self.psiN_RZ(psiv);
        Z0 = self.get('zmaxis');

        #move to the0

        dphi = the0*qv/float(ndphi)
        theta_rot = 0.;

        #initialization of ode solver
        r = ode(self.dr).set_integrator('vode',method = 'adams',with_jacobian = False, 
                max_step = 100000);
        r.set_initial_value([R0,Z0],0.);

        for i in range(ndphi):
            r.integrate(r.t-dphi);

        dphi = 2.*np.pi/ndphi;
        yi = np.array([r.y[0],r.y[1],r.t,the0]);
        for i in range(ndphi*int(qv)):
            r.integrate(r.t-dphi);
            yi = np.vstack((yi,np.array([r.y[0],r.y[1],r.t,the0+float(i+1)*dphi/qv])));
        return yi;


    def transit_path(self, E, mass, charge, sigma,Lambda, psiN, Is_Toroidal=False,
            RZpos = None, N_Max=100000,DTcoef=1., nSig=2,quiet = False):
        """
        Particle orbit during 1 poloidal return
        E is energy in eV unit
        mass is proton ratio
        charge is unit e ratio, proton is 1, electron is -1, alpha is 2
        sigma is parallel or antiparallel
        Lambda is v||^2/v^2
        psi is 0 to 1

        """

        if(RZpos == None): 
            R0, Z0   = self.psin_RZ(psiN,sigma=sigma); 
            phi0 = 0.;
            v    = pb.ev2vi(E,mass);
            vpa  = v*np.sqrt(Lambda)*sigma;
            vpe  = np.sqrt(v**2-vpa**2);
            if not quiet: 
                   print("R0 and Z0 is %f %f"%(R0,Z0));
        else:
            R0 	= RZpos[0];
            phi0= RZpos[1];
            Z0 	= RZpos[2];
            vpa = RZpos[3];
            vpe = RZpos[4];
            v 	= np.sqrt(vpa**2 + vpe**2);
            E   = 0.5*pb.mp*mass*v**2/pb.eV2J;

        Zhe  = charge*pb.eV2J;
        mq   = pb.mp*mass/Zhe;
        mh   = pb.mp*mass;
        B0   = np.sqrt( self.B2(R0,Z0) );
        Rm   = self.data['rmaxis'][0]
        Zm   = self.data['zmaxis'][0]

        mum  = mass*pb.mp*vpe**2/2./B0;
        #dt   = self.data['rmaxis'][0]*2.*np.pi/abs(vpa)/1000.;
        dt   = mq/(B0/np.pi/2.)*DTcoef;
        X0   = np.array([R0,Z0])
        Xt	= np.array([R0,Z0])
        Xm	= np.array([Rm,Zm])
        theN = np.arctan2(np.cross(X0-Xm,Xt-Xm),np.dot(X0-Xm,Xt-Xm))
        if theN<0.: 
              theN += np.pi*2.;

        #calculation for drift frequency
        dXdt = self.dXdt(0,[R0,0.,Z0,vpa,mum],[Zhe,mq,mh]);
        B = self.B_field(R0,Z0);
        Babs = np.sqrt(B[0]**2+B[1]**2+B[2]**2);
        b = B/Babs;
        Vpa = vpa*b;
        Vd= dXdt[0:3]-Vpa;
        q = self.q_RZ(R0,Z0)
        dzeta = np.array([0.,1./R0,0])
        dtheta= self.Grad_Theta(R0,Z0)
        Omegad = np.dot(Vd,dzeta-q*dtheta);

        r    = ode(self.dXdt).set_integrator('vode',method='adams',with_jacobian=False,
                max_step=1000)
        r.set_initial_value([R0,phi0,Z0,vpa,mum],0.)
        r.set_f_params([Zhe,mq,mh]);

        Inds = 0;
        yi   = [[R0,phi0,Z0,vpa]]
        Omeg = [Omegad]
        Sig = 0;

        SS = 0.;
        iind	= 0;
        T_Theta = 0.;
        while Sig<nSig and iind <= N_Max:
            #time evolution
            r.integrate(r.t+dt);

            #calculation for drift frequency
            dXdt = self.dXdt(r.t,r.y,[Zhe,mq,mh]);
            B = self.B_field(r.y[0],r.y[2]);
            Babs = np.sqrt(B[0]**2+B[1]**2+B[2]**2);
            b = B/Babs;
            Vpa = r.y[3]*b;
            Vd= dXdt[0:3]-Vpa;
            q = self.q_RZ(r.y[0],r.y[2]);
            dzeta = np.array([0.,1./r.y[0],0])
            dtheta= self.Grad_Theta(r.y[0],r.y[2])
            Omegad = np.dot(Vd,dzeta-q*dtheta);

            yi += [r.y[0:4]];
            Omeg += [Omegad];

            #End condition check
            theO	= theN
            Xt	= np.array([r.y[0],r.y[2]]);
            theN = np.arctan2(np.cross(X0-Xm,Xt-Xm),np.dot(X0-Xm,Xt-Xm))
            if theN<0:
                 theN+= np.pi*2.
            T_Theta += np.abs(theO - theN);
            if np.abs(theO-theN)  >np.pi and T_Theta > np.pi*2.: 
                 if not quiet:
                     print(theO, theN)
                 Sig += 2;
            #if T_Theta > np.pi*2.:
            #     print(theO, theN)
            #     Sig += 2;
                
                

        yi = np.array(yi)
        Omeg = np.array(Omeg);

        dR = yi[-1,0]-yi[-2,0];
        dp = yi[-1,1]-yi[-2,1];
        dZ = yi[-1,2]-yi[-2,2];
        dV = yi[-2,3]-yi[-2,3];
        dr = (Z0-yi[-2,2])/dZ;
        
        Re = yi[-2,0]+dR*dr;
        pe = yi[-2,1]+dp*dr;
        Ve = yi[-2,3]+dV*dr;
        t = r.t-dt+dt*dr
        yi[-1,0] = Re;
        yi[-1,1] = pe;
        yi[-1,2] = Z0;
        yi[-1,3] = Ve;

        return np.array(yi),t,dt,Omeg;

    def transit_path_RE(self, E, sigma,Lambda,psiN, Is_Toroidal=False, N_Max = 10000000):
        """
        Particle orbit of runawya electrons
        E is energy in eV unit
        Lambda is v||^2/v^2
        sigma is parallel or antiparallel
        psi is 0 to 1

        """

        R0   = self.psiN_RZ(psiN);
        Z0   = self.data['zmaxis'][0];
        phi0 = 0.;
        Zhe  = -pb.eV2J;
        mq   = pb.me/Zhe;
        mh   = pb.me;
        B0   = np.sqrt( self.B2(R0,Z0) );
        Rm   = self.data['rmaxis'][0]
        Zm   = self.data['zmaxis'][0]

        v    = pb.ev2ve(E);
        gamma = 1./np.sqrt(1.-v**2/pb.cv**2);
        vpa  = v*np.sqrt(Lambda)*gamma*sigma;
        vpe  = v*np.sqrt(1-Lambda)*gamma;
        mum  = mh*vpe**2/2./B0;
        #dt   = self.data['rmaxis'][0]*2.*np.pi/abs(vpa)/100.;
        dt   = 1./np.abs(pb.fce_re(B0,E))/100.e0
        theN = np.arctan2(Z0-Zm,R0-Rm);

        mc2 = pb.mec2*pb.eV2J

        #calculation for drift frequency
        dXdt = self.dXdt_RE(0,[R0,0.,Z0,vpa,mum],[Zhe,mq,mh,mc2]);
        B = self.B_field(R0,Z0);
        Babs = np.sqrt(B[0]**2+B[1]**2+B[2]**2);
        b = B/Babs;
        Vpa = vpa*b;
        Vd= dXdt[0:3]-Vpa;
        q = self.q_RZ(R0,Z0)
        dzeta = np.array([0.,1./R0,0])
        dtheta= self.Grad_Theta(R0,Z0)
        Omegad = np.dot(Vd,dzeta-q*dtheta);

        r    = ode(self.dXdt_RE).set_integrator('vode',method='adams',with_jacobian=False,
                max_step=1000)
        r.set_initial_value([R0,phi0,Z0,vpa,mum],0.)
        r.set_f_params([Zhe,mq,mh,mc2]);

        Inds = 0;
        yi   = [[R0,0.,Z0,vpa]]
        Omeg = [Omegad]
        Sig = 0;
        iind = 0;
        #dXdts = [];
        #edrift=[]
        while Sig<2 and iind<=N_Max:
            #time evolution
            r.integrate(r.t+dt);

            #calculation for drift frequency
            dXdt = self.dXdt_RE(r.t,r.y,[Zhe,mq,mh,mc2]);
            #dXdts += [dXdt[0:4]]
            B = self.B_field(r.y[0],r.y[2]);
            Babs = np.sqrt(B[0]**2+B[1]**2+B[2]**2);
            b = B/Babs;
            Vpa = r.y[3]*b;
            Vd= dXdt[0:3]-Vpa;
            q = self.q_RZ(r.y[0],r.y[2]);
            dzeta = np.array([0.,1./r.y[0],0])
            dtheta= self.Grad_Theta(r.y[0],r.y[2])
            Omegad = np.dot(Vd,dzeta-q*dtheta);
            #edrift += [np.dot(self.norm_vec(r.y[0],r.y[2]),self.e_drift(r.y[0],r.y[2],E,Lambda))];
            #vrot += [np.dot(self.norm_vec(r.y[0],r.y[2]),Vpa[:])]
            #End condition check
            if not Is_Toroidal:  
                theO = theN;
                theN = np.arctan2(r.y[2]-Zm,r.y[0]-Rm);
                Inds += 1;
                if theO*theN <0:
                    Sig+=1;
            else:
                if np.abs(r.y[1]) >= np.pi*2.:
                    Sig=2
            yi += [r.y[0:4]];
            Omeg += [Omegad];
            iind +=1

        yi = np.array(yi)
        Omeg = np.array(Omeg);
        

        dR = yi[-1,0]-yi[-2,0];
        dp = yi[-1,1]-yi[-2,1];
        dZ = yi[-1,2]-yi[-2,2];
        dV = yi[-2,3]-yi[-2,3];
        dr = (Zm-yi[-2,2])/dZ;
        
        Re = yi[-2,0]+dR*dr;
        pe = yi[-2,1]+dp*dr;
        Ve = yi[-2,3]+dV*dr;
        t = r.t-dt+dt*dr
        yi[-1,0] = Re;
        yi[-1,1] = pe;
        yi[-1,2] = Zm;
        yi[-1,3] = Ve;

        return np.array(yi),t,dt,Omeg#,np.array(dXdts);

    def transit_Lorbit(self, E, mass, charge, sigma,Lambda, psiN, Is_Toroidal=False,RZpos = None,
            MaxT = None):
        """
        Particle path of Lorentz orbit
        E is energy in eV unit
        mass is proton ratio
        charge is unit e ratio, proton is 1, electron is -1, alpha is 2
        sigma is parallel or antiparallel
        Lambda is v||^2/v^2
        psi is 0 to 1
        MaxT is maximum T of dt unit;


        here dt is 0.01/f where f is cyclotron frequency
        """

        if(RZpos == None): 
            R0   = self.psiN_RZ(psiN); 
            Z0   = self.data['zmaxis'][0]; 
            phi0 = 0.;
            v    = pb.ev2vi(E,mass);
            vpa  = v*np.sqrt(Lambda)*sigma;
            vpe  = np.sqrt(v**2-vpa**2);
        else:
            R0 = RZpos[0];
            phi0 = RZpos[1];
            Z0 = RZpos[2];
            vpa = RZpos[3];
            vpe = RZpos[4];
            v = np.sqrt(vpa**2 + vpe**2);
        if (MaxT == None):
            MaxT = 100000;

        Zhe  = charge*pb.eV2J;
        mq   = pb.mp*mass/Zhe;
        qm   = 1./mq;
        B0   = np.sqrt( self.B2(R0,Z0) );
        omeg = qm*B0;

        rho  = vpe/omeg

        dt   = 1./(omeg/np.pi/2.)*0.01;

        #calculation for drift frequency
        Xyz  = pb.RpZ2xyz(np.array([R0,phi0,Z0]));
        Bxyz = self.Bxyz(Xyz[0],Xyz[1],Xyz[2])
        Xc   = pb.gyration_position(Bxyz,Xyz, vpe, rho, np.pi);
        Vc   = pb.gyration_velocity(Bxyz,Xyz, vpe, rho, np.pi);
        Ri   = pb.xyz2RpZ(Xc);
        Vpe  = pb.Vxyz2RpZ(Xc,Vc);

        Bvec = self.B_field(R0,Z0,cyclic=True);
        Babs = np.sqrt(np.dot(Bvec,Bvec));
        bvec = Bvec/Babs;
        Vpa  = vpa*bvec

        Vi   = Vpa+Vpe;
        Vi[1] /= R0;


        r    = ode(self.dXdt_Lorbit).set_integrator('vode',method='adams',
                with_jacobian=False, max_step=1000)
        r.set_initial_value([Ri[0],Ri[1],Ri[2],Vi[0],Vi[1],Vi[2]],0.)
        r.set_f_params([qm]);

        Inds = 0;
        yi   = [[Ri[0],Ri[1],Ri[2],Vi[0],Vi[1],Vi[2]]];

        while r.t<MaxT*dt:
            #time evolution
            r.integrate(r.t+dt);

            #calculation for drift frequency
            dXdt = self.dXdt_Lorbit(r.t,r.y,[qm]);

            yi += [r.y[0:6]];

        yi = np.array(yi)

        return np.array(yi)

    def transit_Lorbit_RE(self, E, sigma,Lambda, psiN, Is_Toroidal=False,RZpos = None,
            MaxT = None, Init_option = 1):
        """
        Particle path of Lorentz orbit
        E is energy in eV unit
        sigma is parallel or antiparallel
        Lambda is v||^2/v^2
        psi is 0 to 1
        MaxT is maximum T of dt unit;


        here dt is 0.01/f where f is cyclotron frequency
        """

        if(RZpos == None): 
            R0   = self.psiN_RZ(psiN); 
            Z0   = self.data['zmaxis'][0]; 
            phi0 = 0.;
            v    = pb.ev2ve(E);
            vpa  = v*np.sqrt(Lambda)*sigma;
            vpe  = np.sqrt(v**2-vpa**2);
        else:
            R0 = RZpos[0];
            phi0 = RZpos[1];
            Z0 = RZpos[2];
            vpa = RZpos[3];
            vpe = RZpos[4];
            v = np.sqrt(vpa**2 + vpe**2);
        if (MaxT == None):
            MaxT = 100000;

        gamma   = pb.ev2gamma(E);
        Zhe  = -pb.eV2J;
        mq   = pb.me/Zhe;
        qm   = 1./mq;
        B0   = np.sqrt( self.B2(R0,Z0) );
        omeg = qm*B0/gamma;

        rho  = np.abs(vpe / omeg);
        print(rho)

        dt   = 1./np.abs(pb.fce_re(B0, E))*0.01
        #dt   = 1./np.abs(pb.fce(B0))*0.01

        #calculation for drift frequency
        Xyz  = pb.RpZ2xyz(np.array([R0,phi0,Z0]));
        Bxyz = self.Bxyz(Xyz[0],Xyz[1],Xyz[2])
        G_phase = np.pi/2.
        if(Init_option == 1): 
            Xc   = pb.gyration_position(Bxyz,Xyz, vpe, rho, G_phase);
            Vc   = pb.gyration_velocity(Bxyz,Xyz, vpe, rho, G_phase);
            Xpe  = Xc - Xyz; 
        elif(Init_option==2):
            Axyz = np.array([1.,0.,1.]);
            Xc  = np.cross(Bxyz,Axyz);
            Xpe = Xc / np.sqrt(np.dot(Xc,Xc))*rho;
            Xc  = Xyz + Xpe;
            Vc   = np.cross(Xpe,Bxyz)/rho/np.sqrt(np.dot(Bxyz,Bxyz))*vpe

        Bvec = self.B_field(R0,Z0,cyclic=True);
        Babs = np.sqrt(np.dot(Bvec,Bvec));
        icylinder = False;
        if(icylinder):
            Ri   = pb.xyz2RpZ(Xc);
            Vpe  = pb.Vxyz2RpZ(Xc,Vc);
            Vpa  = vpa*Bvec/Babs;
            Vpe[1]  /= Ri[0];
            Vpa[1]  /= Ri[0];

        else: 
            Bxyz = self.Bxyz(Xyz[0],Xyz[1],Xyz[2])
            Vpa = vpa*Bxyz/np.sqrt(np.dot(Bxyz,Bxyz));
            Ri  = Xc;
            Vpe = Vc;

        Vi  = Vpa-Vpe;   # This is for electron rotation direction
        Ui  = Vi*gamma;
        Xi  = Ri;

        #Ui  = Vpa*gamma
        #Xi  = Xyz

        r    = ode(self.dXdt_Lorbit_RE).set_integrator('vode',method='adams',
                with_jacobian=False, max_step=1000)
        r.set_initial_value([Xi[0],Xi[1],Xi[2],Ui[0],Ui[1],Ui[2]],0.)
        r.set_f_params([qm,gamma]);

        Inds = 0;
        yi   = [pb.xyz2RpZ(Xc)[0:3]];
        #yi   = [Xyz[0:3]];

        while r.t<MaxT*dt:
            #time evolution
            r.integrate(r.t+dt);

            #calculation for drift frequency
            dXdt = self.dXdt_Lorbit_RE(r.t,r.y,[qm,gamma]);

            yi += [pb.xyz2RpZ(r.y[0:3])[0:3]];

        yi = np.array(yi)

        return np.array(yi)

    def Init_for_transit_path(self, E, mass, charge, sigma,Lambda, psiN ):
        """
	Initial condition output for transit_path not Lorbit
	Input:
		E	: energy
		mass	: mass in proton mass unit
		charge	: unit in e
		sigma	: parallel and anti parallel
		Lambda	: v||/v
		psiN	: initial psi position from 0 to 1;

        
        """

        R0   = self.psiN_RZ(psiN); 
        Z0   = self.data['zmaxis'][0]; 
        phi0 = 0.;
        v    = pb.ev2vi(E,mass);
        vpa  = v*np.sqrt(Lambda)*sigma;
        vpe  = np.sqrt(v**2-vpa**2);

        return R0, phi0, Z0, vpa, vpe;


    def init_rbx(self): 
       ix1 = 0;
       nbbbs = self.data['nbbbs'][0];
       zmaxis = self.data['zmaxis'][0];
       zbbbs = self.data['zbbbs'][0]
       rbbbs = self.data['rbbbs'][0]
       izblr = np.zeros(2,dtype='int');
       for i in range(nbbbs-1): 
       	if((zmaxis - zbbbs[i])*(zmaxis- zbbbs[i+1])<=0.):
       		izblr[ix1] = i;
       		ix1+=1;
       ####### get rbxl and rblr #######
       #double R2, R1, Z2, Z1, Rb1, Rb2, Zx;
       print("ix1 is %d"%ix1);
       if(ix1==2):
       	Zx  = zmaxis;
       	ix2 = izblr[0];
       	ix3 = izblr[1];
       	R2 = rbbbs[ix2+1]; R1 = rbbbs[ix2];
       	Z2 = zbbbs[ix2+1]; Z1 = zbbbs[ix2];
       	Rb1 = (R2 - R1)/(Z2-Z1)*(Zx-(R2*Z1-R1*Z2)/(R2-R1));
       	R2 = rbbbs[ix3+1]; R1 = rbbbs[ix3];
       	Z2 = zbbbs[ix3+1]; Z1 = zbbbs[ix3];
       	Rb2 = (R2 - R1)/(Z2-Z1)*(Zx-(R2*Z1-R1*Z2)/(R2-R1));
       	print("Rb1 is %f and Rb2 is %f\n"%(Rb1, Rb2));
       	if(Rb1 >= Rb2):
       		self.rbxl = Rb2;
       		self.rbxr = Rb1;
       	else:
       		self.rbxl = Rb1;
       		self.rbxr = Rb2;
       else:
       	print("plasma boundary has problems\n");
       	self.rbxl = rbbbs.min(); 
       	self.rbxr = rbbbs.max();
       
       	##high-field side initialization psi_R_l_spl;
       ddR = 0.02
       nR = 100
       
       rmaxis	= self.data['rmaxis'][0];
       simag	= self.data['simag'][0];
       psiw	= self.data['psiw'][0];
       
       dR = -(rmaxis-self.rbxl+ddR)/float(nR-1);
       Rs = np.zeros(nR,dtype='float');
       psi_Rs = np.zeros(nR,dtype='float');
       for i in range(nR):
       	Rs[i] = rmaxis+dR*float(i);
       	psi_Rs[i] = (self.f(Rs[i],zmaxis)-simag)/psiw;
       
       psi_Rs[0]	= 0.;
       self.f_RZ_Nl   = interpolate.interp1d(psi_Rs,Rs,kind='cubic')
       
       ###low-field side initialization psi_R_l_spl;
       dR = -(rmaxis-self.rbxr-ddR)/float(nR-1);
       for i in range(nR):
       	Rs[i] = rmaxis+dR*float(i);
       	psi_Rs[i] = (self.f(Rs[i],zmaxis)-simag)/psiw;
       psi_Rs[0] 	= 0.; 
       self.f_RZ_Nr   = interpolate.interp1d(psi_Rs,Rs,kind='cubic')
### End of class
	 
if __name__ == '__main__':
    #from pylab import *
    polyg_test = True
    psi_test = True
    q_test = True
    seo_test = True
    fpath = "/home/trhee/equilibrium/21523/"
    filename = fpath + "g021523.005000"

    g1 = geqdsk_dk(filename = filename,gR0B0=False)
    g1.init_operators()
#    g1.get_B_abs()
    nnew = 100
    Rnew = np.linspace(1.4,2.3,100)
    Znew = 0.5
    Bs   = np.zeros(nnew,dtype='float')
    Bvec   = np.zeros([3,nnew],dtype='float')
    b   = np.zeros([3,nnew],dtype='float')
    b2  = np.zeros(nnew,dtype='float')
    b_dot_curl_b   = np.zeros(nnew,dtype='float')
    Curl_b = np.zeros([3,nnew],dtype='float')
    Grad_B = np.zeros([3,nnew],dtype='float')
    for i in range(nnew):
        data1   =g1.EquiB_package_for_dk(Rnew[i],Znew)
        Curl_b[:,i] = data1['curl_b']
        Grad_B[:,i] = data1['Grad_B']

    plt.figure(4)
    plt.plot(Rnew,Curl_b[0,:],'rx-')
    plt.plot(Rnew,Curl_b[1,:],'bx-')
    plt.plot(Rnew,Curl_b[2,:],'kx-')
    plt.legend()
    plt.xlabel('R')
    plt.title('curl_b')


    plt.figure(7)
    plt.plot(Rnew,Grad_B[0,:],'rx')
    plt.plot(Rnew,Grad_B[2,:],'bx')
    plt.legend()
    plt.xlabel('R')


    psi_n = 0.8
    g1 = geqdsk_dk(filename = filename,gR0B0=False)
    g1.init_data_o(g1.BtBphi);  #initialize the surface data

    ndphi = 50


    if psi_test:
        for psi_n in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.992]:
            averv = g1.surf_aver(psi_n,ndphi = ndphi)
            val = g1.get('simag')+g1.get('psiw')*psi_n;
            print("Psi %f Error is %f"%(psi_n, averv/val - 1.))

    if polyg_test:
        ndphi = 50;
        psi_n = 0.8;
        y = g1.polyg(psi_n,ndphi = ndphi);
        leny = y.shape[0]
        for i in range(leny):
            print(y[i,0],y[i,1],y[i,2]/np.pi/2.,np.abs(g1.BtBphi(y[i,0],y[i,1])))
        plt.figure(1);
        plt.plot(y[:,0],y[:,1],'ro-');
        plt.plot(y[0,0],y[0,1],'bd-',markersize=10);
        plt.plot(y[-1,0],y[-1,1],'rx-',markersize=10);
        plt.plot(g1.get('rbbbs'),g1.get('zbbbs'),'b-');
        print("q is %f"%(g1.q_inter(psi_n)));
        plt.show();
    #q surface test
    if q_test:
        ndphi = 100
        g1.init_data_o(g1.BtBphi);  #initialize the surface data
        for psi_n in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
            averv = g1.surf_aver_q(psi_n,ndphi = ndphi)
            val = g1.q_inter(psi_n);
            print("Psi %f Calculated %f equilibrium  %f error %f"%(psi_n, averv, val, val/averv-1.));
    if seo_test:
        for qv in [4, 4.5, 5, 5.5, 6]:
            psiv = g1.q_root_inter(qv)
            R,Z = g1.psin_RZ(psiv);
            psit = g1.f_normal(R,Z)[0]
            print("At q=%f error of calculated psi value is %f"%(qv, (psiv-psit)/psiv*100.))
        yi = g1.Closed_Lines(5, np.pi/2., ndphi = 50)

        plt.figure(1)
        plt.plot(yi[:,0],yi[:,1],'ro-')

        fig1 = plt.figure(2)
        ax = fig1.add_subplot(111,projection='3d');
        ax.plot(yi[:,0]*np.cos(yi[:,2]),yi[:,0]*np.sin(yi[:,2]),yi[:,1])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        print(yi[:,3]/np.pi/2.)

        plt.show()

    plt.show()


