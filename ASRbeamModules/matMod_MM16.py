# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:07:30 2019

@author: simkon
"""

import numpy as np

def matCalc_concrete(eps0, deps, time, dtime, predef, dpredef, matProps, sig, stiff, stateVars):
    
    E_0, ft, phi, fc, Gc, alpha_critical, sig_u, sig_L, beta_E = matProps[:9] 
    
    sig, stiff, stateVars = np.copy(sig), np.copy(stiff), np.copy(stateVars)
    
    eps_asr = stateVars[0]
    deps_asr =   depsAsrCalc(sig, dpredef, sig_u, sig_L) 
    eps_asr = eps_asr + deps_asr
    
    
    
    eps = eps0 + deps
    
    eps_el = eps - eps_asr
    E = E_0/(1.+phi)
    sig = E * eps_el
    
    stiff = E
    
    stateVars[0] = eps_asr
    
        
    return sig, stiff, stateVars

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
        W = (sig-sig_u)/(sig_L-sig_u)   
    else:
        W = 1.
    if W < 0.:
        print('ERROR: NEGATIVE W')

    deps_asr  = deps_asr_free * W
    
    return deps_asr
