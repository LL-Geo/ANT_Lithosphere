import numpy as np
from math import erf
from scipy import optimize

# Input parameters
# Gas constant
R=8.3145
# Reference pressure
Pr=1.5e9
# Reference temperature
TKr=1473.
T0=273.
# Grain size and reference grain size
d=1.e-3
dr=1.e-3
# Reference density
rhor=3176.15
# Thermal expansivity
alphaT=3.59e-5
# Bulk modulus
bmod=115.2
# Raw frequency
freq=0.01

# Brent temperature minimization bounds
AX=0.
CX=2000.
tol=1e-3

# Other parameters (density, compressibility etc.):
rho0=3300.
a0=2.832e-5
a1=0.758e-8
K0=130.e9
KT=4.8
grun=6.
p0=3330.
AY=1.0
CY=1.6

# Set key parameters
Ab=0.664
alpha=0.38
tauP=6.e-5
Teta=0.94
beta=0.
delphi=0.
gamma=5.
lambdaphi=0.

# Set dV0 vs. P parameterisation parameters
y_a = -7.334963115431564676e-23
y_b = 7.510867653681621105e-12
y_c = 1.000184023114681908e+00

# Set solidus temperature at 50 km depth
solidus_50km = np.loadtxt('./data/potential_temperature/solidus_50km_temperature.T').astype(float)

def funcVs(T,Vs_obs,m,dep,Mg):

  # Difference between observed and calculated temperature
  funcVs=np.abs(Vs_obs-Vs_calc(m,T,dep,Mg))
  return funcVs

def funcV0(x,Pin,K0,KT):

  # Difference between observed and calculated pressure for given value of P
  funcV0=np.abs(((K0*(3./2.)*(x**(7./3.)-x**(5./3.))*(1.+(((3./4.)*(KT-4.))*(x**(2./3.)-1.))))*1e-9)-(Pin*1e-9))
  return funcV0

def Vs_calc(m,T,dep,Mg):

  mu0 = m[0]
  dmudT = m[1]
  dmudP = m[2]
  eta0 = 10**m[3]
  E = m[4]
  Va = m[5]
  dTdz = m[6]
  dmudMg = m[7]
  dpdMg = m[8]
  sol50 = solidus_50km

  Pg=(dep/30.)
  P=Pg*1.e9
  TK=T+273.

  Tsol=sol50+(dTdz*(dep-50.))
  Tn=TK/(Tsol+273.)
  # Initialise parameters for raw Vs
  if Tn < Teta:
    Aeta=1.
  elif Tn >= Teta and Tn<1.:
    Aeta=np.exp((-1.*((Tn-Teta)/(Tn-(Tn*Teta))))*np.log(gamma))
  else:
    Aeta=(1./gamma)*np.exp(-delphi)

  # Work out viscosity given A
  eta=((eta0*np.exp((E/R)*(1./TK-1./TKr))*np.exp((Va/R)*(P/TK-Pr/TKr)))*Aeta)

  # Unrelaxed compliance
  Ju=1./(1.e9*(mu0+(dmudP*Pg)+(dmudT*T)+(dmudMg*(100-Mg)/100)))

  # Determine input parameters for complex compliance terms
  tauM=eta*Ju
  tau=(3.*dep)/4.2
  tauS=tau/(2*np.pi*tauM)
  if Tn < 0.91:
    Ap=0.01
  elif Tn>=0.91 and Tn<0.96:
    Ap=0.01+(0.4*(Tn-0.91))
  elif Tn>=0.96 and Tn<1.:
    Ap=0.03
  else:
    Ap=0.03+beta
  if Tn < 0.92:
    sigmap=4.
  elif Tn>=0.92 and Tn<1.:
    sigmap=4.+(37.5*(Tn-0.92))
  else:
    sigmap=7.

  # Determine complex compliance terms
  J1=Ju*(1.+((Ab*(tauS**alpha))/alpha)+((np.sqrt(2.*np.pi)/2.)*Ap*sigmap*(1.-erf((np.log(tauP/tauS))/(np.sqrt(2.)*sigmap)))))

  # include pressure and temperature-dependent alpha
  #dV0=optimize.brent(funcV0,brack=(AY,CY),args=(P,K0,KT,),tol=tol)
  dV0=y_a*P**2 + y_b*P + y_c
  alphaP0=dV0*np.exp((grun+1.)*((dV0**(-1.))-1.))
  rhoP0=(p0+dpdMg*(Mg-89.27)/100)*dV0
  intalphaT=(a0*(TK-273.))+((a1/2.)*((TK**2.)-(273.**2.)))
  rho=rhoP0*(1.-(alphaP0*intalphaT))

  # Calculate Vs
  Vs=1./(np.sqrt(rho*J1)*1000.)
  return Vs

def T_calc(m,Vs,dep,Mg):

  # Calculate temperature from Vs based on optimisation
  T=optimize.brent(funcVs,brack=(AX,CX),args=(Vs,m,dep,Mg,),tol=1e-3)

  return T

def Q_calc(m,Vs,dep,Mg):

  mu0 = m[0]
  dmudT = m[1]
  dmudP = m[2]
  eta0 = 10**m[3]
  E = m[4]
  Va = m[5]
  dTdz = m[6]
  sol50 = solidus_50km
  dmudMg = m[7]
  dpdMg = m[8]
  # Initialise parameters for attenuation
  T=optimize.brent(funcVs,brack=(AX,CX),args=(Vs,m,dep,Mg,),tol=tol)
  TK=T+273.
  Pg=(dep/30.)
  P=Pg*1.e9

  Tsol=sol50+(dTdz*(dep-50.))
  Tn=TK/(Tsol+273.)

  # Initialise parameters for raw Vs
  if Tn < Teta:
    Aeta=1.
  elif Tn >= Teta and Tn<1.:
    Aeta=np.exp((-1.*((Tn-Teta)/(Tn-(Tn*Teta))))*np.log(gamma))
  else:
    Aeta=(1./gamma)*np.exp(-delphi)

  # Work out viscosity given A
  eta=((eta0*np.exp((E/R)*(1./TK-1./TKr))*np.exp((Va/R)*(P/TK-Pr/TKr)))*Aeta)

  # Unrelaxed compliance
  Ju=1./(1.e9*(mu0+(dmudP*Pg)+(dmudT*T)+(dmudMg*(100-Mg)/100)))

  # Determine input parameters for complex compliance terms
  tauM=eta*Ju
  tau=(3.*dep)/4.2
  tauS=tau/(2*np.pi*tauM)
  if Tn < 0.91:
    Ap=0.01
  elif Tn>=0.91 and Tn<0.96:
    Ap=0.01+(0.4*(Tn-0.91))
  elif Tn>=0.96 and Tn<1.:
    Ap=0.03
  else:
    Ap=0.03+beta
  if Tn < 0.92:
    sigmap=4.
  elif Tn>=0.92 and Tn<1.:
    sigmap=4.+(37.5*(Tn-0.92))
  else:
    sigmap=7.

  # Determine complex compliance terms
  J1=Ju*(1.+((Ab*(tauS**alpha))/alpha)+((np.sqrt(2.*np.pi)/2.)*Ap*sigmap*(1.-erf((np.log(tauP/tauS))/(np.sqrt(2.)*sigmap)))))
  J2=(Ju*(np.pi/2.)*((Ab*(tauS**alpha))+(Ap*(np.exp(-1.*(((np.log(tauP/tauS))**2.)/(2.*(sigmap**2.))))))))+(Ju*tauS)

  # Calculate Q
  Q=J2/J1
  return Q

def visc_calc(m,Vs,dep,Mg):

  mu0 = m[0]
  dmudT = m[1]
  dmudP = m[2]
  eta0 = 10**m[3]
  E = m[4]
  Va = m[5]
  dTdz = m[6]
  sol50 = solidus_50km
  dmudMg = m[7]
  dpdMg = m[8]  
  # Initialise parameters for viscosity
  T=optimize.brent(funcVs,brack=(AX,CX),args=(Vs,m,dep,Mg,),tol=tol)
  TK=T+273.
  Pg=(dep/30.)
  P=Pg*1.e9

  # Find viscosity
  Tsol=sol50+(dTdz*(dep-50.))
  Tn=TK/(Tsol+273.)

  # Initialise parameters for raw Vs
  if Tn < Teta:
    Aeta=1.
  elif Tn >= Teta and Tn<1.:
    Aeta=np.exp((-1.*((Tn-Teta)/(Tn-(Tn*Teta))))*np.log(gamma))
  else:
    Aeta=(1./gamma)*np.exp(-delphi)

  # Work out viscosity
  eta=((eta0*np.exp((E/R)*(1./TK-1./TKr))*np.exp((Va/R)*(P/TK-Pr/TKr)))*Aeta)
  return eta
