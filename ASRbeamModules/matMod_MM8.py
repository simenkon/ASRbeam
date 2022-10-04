# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:07:30 2019

@author: simkon
"""

import numpy as np

def matCalc_concrete(eps0, deps, time, dtime, predef, dpredef, matProps, stress, stiff, stateVars):
    
    E_0, ft, gamma, fc, Gc, alpha_critical, sig_u, sig_L, beta_E = matProps[:9] 
    E_creep = matProps[9:14]
    lambda_creep_array = matProps[14:19]
    beta_Ecreep = np.array([[-0.01914656,  0.07135232,  0.26914982,  0.227365  ,  0.11070075],
                            [-0.05485107,  0.11590657,  0.30621064,  0.26957512,  0.27000593],
                            [-0.05014045,  0.1123435 ,  0.24394659,  0.23264344,  0.18283563],
                            [-0.0339299 ,  0.09002191,  0.28048416,  0.23998737,  0.18048002],
                            [-0.0555953 ,  0.11949636,  0.24336659,  0.24786201,  0.15465828]])
    stress, stiff, stateVars = np.copy(stress), np.copy(stiff), np.copy(stateVars)
    alpha, eps_asr = stateVars[:2]
    zeta, xi = stateVars[4:9], stateVars[9:14]

    N_creep = len(E_creep)
    
    
    
    #--STRAIN CALCULATIONS--
    deps_asr = depsAsrCalc(stress, dpredef, sig_u, sig_L) 

    eps_asr = eps_asr + deps_asr
    eps_creep = epsCreepCalc(lambda_creep_array, zeta, xi, dtime)

    eps = eps0 + deps
    eps_i = eps - eps_asr - eps_creep
    
    Ds = (1. - eps_asr/(eps_asr+beta_E))*E_0
    sig = Ds * eps_i

    stiff = Ds
    dsig = sig - stress
    stress  = sig
    stateVars[1] = eps_asr
    
    time_mid = time + dtime/2.
    for k in range(N_creep):
        Ccreep_k = creepCompliance_calc(lambda_creep_array, E_creep[k], beta_Ecreep[k], time_mid)
        zeta[k] = zeta[k] + Ccreep_k *dsig
        xi[k]   = np.exp(-dtime/lambda_creep_array[k]) * xi[k] + lambda_creep_array[k]*Ccreep_k/dtime*(1-np.exp(-dtime/lambda_creep_array[k]))*dsig
    

    stateVars[4:9], stateVars[9:14] = zeta, xi

    return stress, stiff, stateVars

def matCalc_steel(eps0, deps, time, dtime, predef, dpredef, matProps,sig, stiff,stateVars):
    
    E_0, f_y, sh = matProps 
    sig, stiff, stateVars = np.copy(sig), np.copy(stiff), np.copy(stateVars)
    
    eps = eps0 + deps
    eps_pl = 0. #plasticity
    
    eps_el = eps - eps_pl
    
    sig = E_0 * eps_el
    
    stiff = E_0
    
        
    return sig, stiff, stateVars


def depsAsrCalc(sig, dpredef, sig_u, sig_L):
    
    #-----------ASR CALCULATIONS---------------
    
    deps_asr_free = dpredef

    if sig < sig_u:
        W = 0.
    elif sig < sig_L:
        W = 1. - np.log10(sig/sig_L)/np.log10(sig_u/sig_L)
    else:
        W = 1.
    if W < 0.:
        print('ERROR: NEGATIVE W')

    deps_asr  = deps_asr_free * W
    
    return deps_asr


def epsCreepCalc(lambda_array, zeta, xi, dtime):
    N_creep = len(lambda_array)
    eps_creep = 0.
    for k in range(N_creep):
        eps_creep = eps_creep + zeta[k] - np.exp(-dtime/lambda_array[k])*xi[k]
    
    return eps_creep


def creepCompliance_calc(lambda_array, Ecreep_k, beta_Ecreep_k, tau):
    
    N_creep = len(lambda_array)
      
    Ccreep_k = 1.0
    for j in range(N_creep):
        Ccreep_k = Ccreep_k - beta_Ecreep_k[j] * (1. - np.exp(-tau/lambda_array[j]) )

    Ccreep_k = Ccreep_k * 1./Ecreep_k
    return  Ccreep_k
	   
    