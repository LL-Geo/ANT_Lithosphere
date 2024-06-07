import numpy as np
from thermodynamic_conversions import *

def likelihood_Vs_plate(data, m, h):
    # Load in file of structure Vs, sig_Vs, T, depth w/ header 
    Vs = data[0]
    sig_Vs = data[1]
    weighted_sig_Vs = (10**h)*sig_Vs
    T = data[2]
    depth = data[3]
    Mg = data[4]
    n = len(T)
    Vs_fromT = np.zeros(n)
    Vs_diff = np.zeros(n)
    RMS = 0 
    for i in range(n):
        Vs_fromT[i] = Vs_calc(m, T[i], depth[i], Mg[i])
        Vs_diff[i] = ((Vs[i] - Vs_fromT[i]) / weighted_sig_Vs[i])**2
        RMS += (Vs[i] - Vs_fromT[i])**2
    chi_squared = np.sum(Vs_diff)
    P = -0.5 * chi_squared - np.log(((2*np.pi)**(n/2))*np.prod(weighted_sig_Vs, dtype=np.longdouble)*1e+150) + np.log(1e+150)
    RMS = np.sqrt(RMS / n) / np.mean(sig_Vs)
    return P, RMS

def likelihood_Vs_adiabat(data, m, h):
    # Load in file of structure Vs, sig_Vs, T, depth w/ header 
    Vs = data[0]
    sig_Vs = data[1]
    weighted_sig_Vs = (10**h)*sig_Vs
    T = data[2]
    depth = data[3]
    Mg = data[4]    
    n = len(T)
    Vs_fromT = np.zeros(n)
    Vs_diff = np.zeros(n)
    RMS = 0 
    for i in range(n):
        Vs_fromT[i] = Vs_calc(m, T[i], depth[i], Mg[i])
        Vs_diff[i] = ((Vs[i] - Vs_fromT[i]) / weighted_sig_Vs[i])**2
        RMS += (Vs[i] - Vs_fromT[i])**2
    chi_squared = np.sum(Vs_diff)
    P = -0.5 * chi_squared - np.log(((2*np.pi)**(n/2))*np.prod(weighted_sig_Vs, dtype=np.longdouble)*1e+150) + np.log(1e+150)
    RMS = np.sqrt(RMS / n) / np.mean(sig_Vs)
    return P, RMS

def likelihood_Vs_adiabat_t(data, m, h):
    # Load in file of structure Vs, T, sig_T, depth w/ header 
    Vs = data[0]
    T = data[2]
    sig_T = np.full(np.shape(T), 100)
    weighted_sig_T = (10**h)*sig_T
    depth = data[3]
    Mg = data[4]  
    n = len(T)
    T_fromVs = np.zeros(n)
    T_diff = np.zeros(n)
    RMS = 0 
    for i in range(n):
        T_fromVs[i] = T_calc(m, Vs[i], depth[i],Mg[i])
        T_diff[i] = ((T[i] - T_fromVs[i]) / weighted_sig_T[i])**2
        RMS += (T[i] - T_fromVs[i])**2
    chi_squared = np.sum(T_diff)
    P = -0.5 * chi_squared - np.log(((2*np.pi)**(n/2))*np.prod(weighted_sig_T, dtype=np.longdouble))
    RMS = np.sqrt(RMS / n) / np.mean(sig_T)
    return P, RMS

def likelihood_attenuation(data, m, h):
    # Load in file of structure Q, sig_Q, Vs, depth w/ header 
    Q = data[0]
    sig_Q = data[1]
    weighted_sig_Q = (10**h)*sig_Q
    Vs = data[2]
    depth = data[3]
    Mg = data[4]
    n = len(Vs)
    Q_fromVs = np.zeros(n)
    Q_diff = np.zeros(n)
    RMS = 0
    for i in range(n):
        Q_fromVs[i] = Q_calc(m, Vs[i], depth[i], Mg[i])
        Q_diff[i] = ((Q[i] - Q_fromVs[i]) / weighted_sig_Q[i])**2
        RMS += (Q[i] - Q_fromVs[i])**2
    chi_squared = np.sum(Q_diff)
    P = -0.5 * chi_squared - np.log(((2*np.pi)**(n/2))*np.prod(weighted_sig_Q, dtype = np.longdouble))
    RMS = np.sqrt(RMS / n) / np.mean(sig_Q)
    return P, RMS

def likelihood_viscosity(data, m, h):
    # Load in file of structure eta, sig_eta, T, depth w/ header
    eta = data[0]
    sig_eta = data[1]
    weighted_sig_eta = (10**h)*sig_eta
    Vs = data[2]
    depth = data[3]
    Mg = data[4]
    n = len(Vs)
    eta_fromVs = np.zeros(n)
    for i in range(n):
        eta_fromVs[i] = np.log10(visc_calc(m, Vs[i], depth[i], Mg[i]))
    mean_eta = np.mean(eta)
    mean_eta_fromVs = np.mean(eta_fromVs)
    chi_squared = ((mean_eta - mean_eta_fromVs) / weighted_sig_eta[0])**2
    P = -0.5 * chi_squared - np.log(((2*np.pi)**(1/2))*weighted_sig_eta[0])
    RMS = np.sqrt((mean_eta - mean_eta_fromVs)**2) / sig_eta[0]
    return P, RMS


def likelihood(data, m, h, n_plate, n_adiabat, n_attenuation, n_viscosity):
    P_xenolith = 0
    P_plate = 0
    P_adiabat = 0
    P_attenuation = 0
    P_viscosity = 0
    RMS = np.zeros(n_plate + n_adiabat + n_attenuation + n_viscosity)
    if n_plate > 0:
        P_plate, RMS[0] = likelihood_Vs_plate(data[0][0], m, h[0])
    if n_adiabat > 0:   
        P_adiabat, RMS[n_plate] = likelihood_Vs_adiabat_t(data[1][0], m, h[n_plate])
    if n_attenuation > 0:
        P_attenuation, RMS[n_plate + n_adiabat] = likelihood_attenuation(data[2][0], m, h[n_plate + n_adiabat])
    if n_viscosity > 0:
        P_viscosity, RMS[n_plate + n_adiabat + n_attenuation] = likelihood_viscosity(data[3][0], m, 0)
    return P_xenolith + P_plate + P_adiabat + P_attenuation + P_viscosity, RMS