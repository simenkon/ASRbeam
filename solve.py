# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:18:53 2020

@author: simkon
"""


#Main program
if __name__ == '__main__':
    pass
    
import numpy as np
import matplotlib.pyplot as plt

import ASRbeamModules.functions as fu
import ASRbeamModules.matMod_MM1 as mat1
import ASRbeamModules.matMod_MM2 as mat2
import ASRbeamModules.matMod_MM3 as mat3
import ASRbeamModules.matMod_MM4 as mat4
import ASRbeamModules.matMod_MM5 as mat5
import ASRbeamModules.matMod_MM6 as mat6
import ASRbeamModules.matMod_MM7 as mat7
import ASRbeamModules.matMod_MM8 as mat8
import ASRbeamModules.matMod_MM9 as mat9
import ASRbeamModules.matMod_MM10 as mat10
import ASRbeamModules.matMod_MM11 as mat11
import ASRbeamModules.matMod_MM12 as mat12
import ASRbeamModules.matMod_MM13 as mat13
import ASRbeamModules.matMod_MM14 as mat14
import ASRbeamModules.matMod_MM15 as mat15
import ASRbeamModules.matMod_MM16 as mat16
import ASRbeamModules.matMod_MM17 as mat17
import ASRbeamModules.matMod_MM18 as mat18
import ASRbeamModules.matMod_MM19 as mat19
import ASRbeamModules.matMod_MM20 as mat20



import ASRbeamModules.preProcess as pre
import ASRbeamModules.postProcess as post

import argparse
import datetime

# Create the parser
parser = argparse.ArgumentParser()
# Add an argument
parser.add_argument('--aname', type=str, required=False)
parser.add_argument('--mat', type=int, required=False)
# Parse the argument
args = parser.parse_args()

# Assign
if args.aname:
    analysis_name = args.aname
else:
    current_time = datetime.datetime.now()
    analysis_name = f'A_{current_time.year}_{current_time.month}_{current_time.day}_{current_time.hour}_{current_time.minute}'

if args.mat:
    MM = args.mat
else:
    print("""No material id number was specified by --mat,
so the default material id = 1 was used, i.e. ASRbeamModules.matMod_MM1""")
    MM = 1


if MM == 1:
    mat = mat1
elif MM == 2:
    mat = mat2
elif MM == 3:
    mat = mat3
elif MM == 4:
    mat = mat4
elif MM == 5:
    mat = mat5
elif MM == 6:
    mat = mat6
elif MM == 7:
    mat = mat7
elif MM == 8:
    mat = mat8
elif MM == 9:
    mat = mat9
elif MM == 10:
    mat = mat10
elif MM == 11:
    mat = mat11
elif MM == 12:
    mat = mat12
elif MM == 13:
    mat = mat13
elif MM == 14:
    mat = mat14
elif MM == 15:
    mat = mat15
elif MM == 16:
    mat = mat16
elif MM == 17:
    mat = mat17
elif MM == 18:
    mat = mat18
elif MM == 19:
    mat = mat19
elif MM == 20:
    mat = mat20
else:
    print('ERROR: Could not find the material model.')
    exit
    


#-PREPROCESSING
#--READING FILES
nodes, elements, crossSections, mater_c, mater_s, loads_step1, loads_step2, analysis, predef_field_1, predef_field_2 = pre.readin('INPUT/')
#--END READING FILES

#--CREATE ARRAYS:
#---CONNECTIVITY MATRIX
C = pre.connectivityMat(nodes, elements)
#---BOUNDARY CONDITIONS
BC = pre.boundaryCondMat(nodes)
#---ELEMENT COORDINATES
elCoord = pre.elCoordMat(nodes, elements)
#---CROSS SECTION ARRAY
crossSection_array = pre.crossSectionMat(crossSections)
#---MATERIAL PROPERTIES
#---FOR NOW I HAVE ONLY INCLUDED ONE MATERIAL FOR CONCRETE AND ONE FOR STEEL  
NmatProps_c = int(mater_c[0][1])
NstateVars_c = int(mater_c[0][2])
NmatProps_s = int(mater_s[0][1])
NstateVars_s = int(mater_s[0][2])
matProps_c = mater_c[0][3:3+NmatProps_c]
matProps_s = mater_s[0][3:3+NmatProps_s]

#--NUMBER OF ELEMENTS AND DOFS
Nel = np.shape(elements)[0]
Nnodes = np.shape(nodes)[0]
Ndof = 3*Nnodes + Nel

#-END PREPROCESSING




#--ANALYSIS PARAMETERS
totTime_step1 = analysis[0,0]
Ninc_step1 = int(analysis[0,1])
totTime_step2 = analysis[1,0]
Ninc_step2 = int(analysis[1,1])
Ninc = Ninc_step1 + Ninc_step2 

time_array_step1 = np.linspace(0.,totTime_step1, Ninc_step1+1)
time_array_step2 = np.linspace(0.,totTime_step2, Ninc_step2+1) # np.logspace(np.log(1.),np.log(totTime_step2+1),Ninc_step2+1, base = np.e)-1. #
time_array = np.append(time_array_step1, time_array_step1[-1] + time_array_step2[1:])

#---MAXIMUM NUMBER OF ITERATIONS FOR NEWTON RAPHSON
Nitr = 100


#-HISTORY VARIABLES
#--DISPLACEMENT
D = np.zeros((Ninc+1, Ndof), dtype=np.float64) 

#--CONCRETE STRESS, STRAIN AND STATE VARIABLES
sig_c = np.empty(Nel, dtype=object)
eps_c = np.empty(Nel, dtype=object)
stateVars_c = np.empty(Nel, dtype=object)

for e in range(Nel):
    crossSection_flag = elements[e][1]
    Nc =  len(crossSection_array[crossSection_flag]['zc'])
    sig_c[e] = np.zeros((Ninc+1,2,Nc), dtype=np.float64)
    eps_c[e] = np.zeros((Ninc+1,2,Nc), dtype=np.float64)
    stateVars_c[e] = np.zeros((Ninc+1,2,Nc,NstateVars_c), dtype=np.float64)


#--SET TS_bool FOR EACH FIBRE
for e in range(Nel):
    crossSection_flag = elements[e][1]
    Nc =  len(crossSection_array[crossSection_flag]['zc'])
    crossSection_zc = crossSection_array[crossSection_flag]['zc']    
    for inc in range(Ninc+1):
        for i in range(2):
            for j in range(Nc):
                if crossSection_zc[j] < crossSections[crossSection_flag,6]: 
                    stateVars_c[e][inc][i][j][3] = False
                else:
                    stateVars_c[e][inc][i][j][3] = True


#--STEEL STRESS, STRAIN AND STATE VARIABLES
sig_s = np.empty(Nel, dtype=object)
eps_s = np.empty(Nel, dtype=object)
stateVars_s = np.empty(Nel, dtype=object)
for e in range(Nel):
    crossSection_flag = elements[e][1]
    Ns_layers =  int(crossSections[crossSection_flag,3])
    sig_s[e] = np.zeros((Ninc+1,2,Ns_layers), dtype=np.float64)
    eps_s[e] = np.zeros((Ninc+1,2,Ns_layers), dtype=np.float64)
    stateVars_s[e] = np.zeros((Ninc+1,2,Ns_layers,NstateVars_s), dtype=np.float64)
    


#-EXTERNAL LOAD VECTOR ARRAY    
Rext_tot_step1 = np.zeros(Ndof, dtype=np.float64)
Rext_tot_step2 = np.zeros(Ndof, dtype=np.float64)
for e in range(Nel):
    Rext_el = fu.elExternalForce(elCoord[e],loads_step1[e][1:3])
    for i in range(7):
        Rext_tot_step1[C[e][i]] = Rext_tot_step1[C[e][i]] + Rext_el[i]
    Rext_el = fu.elExternalForce(elCoord[e],loads_step2[e][1:3])
    for i in range(7):
        Rext_tot_step2[C[e][i]] = Rext_tot_step2[C[e][i]] + Rext_el[i]

        
Rext_array = np.zeros((Ninc+1, Ndof), dtype=np.float64)
load_prop_array_step1 = time_array_step1/time_array_step1[-1]
for i in range(1,Ninc_step1+1):
    Rext_array[i] = load_prop_array_step1[i]*Rext_tot_step1
load_prop_array_step2 = time_array_step2/time_array_step2[-1]
for i in range(1, Ninc_step2+1):
    Rext_array[Ninc_step1+i] = Rext_tot_step1 + load_prop_array_step2[i]*(Rext_tot_step2 - Rext_tot_step1)

#-PREDEFINED FIELD ARRAY
def S_func(t, tau_c, tau_l):
    s_value = (1-np.exp(-t/tau_c))/(1+np.exp(-t/tau_c+tau_l/tau_c))
    return s_value
predef_array = np.zeros((Ninc+1,Nel,2))
#for i in range(1,Ninc_step1+1):
#    predef_array[i] = predef_field_1[:,1:3]
load_prop_array_step2 = time_array_step2/time_array_step2[-1]
for i in range(1,Ninc_step2+1):
    predef_array[Ninc_step1+i] = load_prop_array_step2[i] * predef_field_2[:,1:3]



#--INITIALIZE THE MATERIAL STIFFNESS MATRICES FOR EACH ELEMENT, GAUSS POINT, AND SIMPSON POINT
St_c = np.ones((Nel, 2, Nc), dtype=float)
St_s = np.ones((Nel, 2, Ns_layers), dtype=float)



#-SOLVE------------------------------------------------------------------------
#--------------------------------
# INCREMENTAL ITERATIVE PROCEDURE FOR A STEP
#--------------------------------

for inc in range(Ninc):
    converegedIncrement = False
    print('INCREMENT {}'.format(inc))


    time = time_array[inc]
    dtime= time_array[inc+1]-time_array[inc]
    
    #-------------------------
    #-EXTERNAL FORCES, Rext (THIS VECTOR IS KNOWN FOR EACH TIME STEP)
    #-------------------------
    Rext = Rext_array[inc+1]
    
    #-------------------------
    #-INTERNAL FORCES, Rint, AT THE END OF THE TIME INCREMENT GIVEN NO CHANGE IN DISPLACEMENT
    #-------------------------
    #--CALL MATERIAL SUBROUTINES TO UPDATE, STIFFNESS MATRIX, AND STATE VARS
    
    deps = 0.
    for e in range(Nel):
        crossSection_flag = elements[e][1]
        Nc_layers =  2*int(crossSections[crossSection_flag,2]) + 1
        crossSection_zc = crossSection_array[crossSection_flag]['zc']
        for i in range(2):
            for j in range(Nc_layers):
                predef = predef_array[inc][e][0] + crossSection_zc[j]*predef_array[inc][e][1]
                dpredef = predef_array[inc+1][e][0] + crossSection_zc[j]*predef_array[inc+1][e][1]-predef
                sig_c[e][inc+1][i][j], St_c[e][i][j], stateVars_c[e][inc+1][i][j] =   mat.matCalc_concrete(eps_c[e][inc][i][j], deps, time, dtime, predef, dpredef, matProps_c, sig_c[e][inc][i][j], St_c[e][i][j], stateVars_c[e][inc][i][j]) # mat_f.matcalc_concrete(eps_c[e][inc][i][j], deps, time, dtime, predef, dpredef, matProps_c, sig_c[e][inc][i][j], stateVars_c[e][inc][i][j])

    for e in range(Nel):
        crossSection_flag = elements[e][1]
        Ns_layers =  int(crossSections[crossSection_flag,3])
        for i in range(2):
            for j in range(Ns_layers):            
                sig_s[e][inc+1][i][j], St_s[e][i][j], stateVars_s[e][inc+1][i][j]  = mat.matCalc_steel(eps_s[e][inc][i][j], deps, time, dtime, predef, dpredef, matProps_s, sig_s[e][inc][i][j], St_s[e][i][j], stateVars_s[e][inc][i][j]) # mat_f.matcalc_steel(eps_s[e][inc][i][j], deps, time, dtime, predef, dpredef, matProps_s, sig_s[e][inc][i][j],   stateVars_s[e][inc][i][j])# #deps = 0.
    
    
    Rint = np.zeros(Ndof, dtype=float)
    for e in range(Nel):
        crossSection_flag = elements[e][1]
        crossSection_el = crossSection_array[crossSection_flag]
        Rint_el = fu.elInternalForce(crossSection_el,elCoord[e], sig_c[e][inc+1], sig_s[e][inc+1])
        for i in range(7):
            Rint[C[e][i]] = Rint[C[e][i]] + Rint_el[i]
    
    
    #-------------------------
    #-OUT OF BALANCE FORCE, R, difference between external forces Rext and internal forces Rint
    #-------------------------
    R = Rext - Rint
    #-------------------------
    
    #-------------------------
    #-ITERATIVE PROCEDURE TO FIND DISPLACEMENT dD FOR THE TIME INCREMENT FROM t_n TO t_(n+1)
    #--The displacement increment dD is initialized to the zero vector, and updated with ddD for every iteration
    #-------------------------
    dD = np.zeros(Ndof, dtype=float)        
    for itr in range(Nitr):
        #start = test_time.time()
        #-UPDATE THE GLOBAL STIFFNESS MATRIX
        K = np.zeros((Ndof,Ndof), dtype=float)
        for e in range(Nel):
            crossSection_flag = elements[e][1]
            crossSection_el = crossSection_array[crossSection_flag]
            K_e = fu.elStiffMat(crossSection_el, elCoord[e], St_c[e], St_s[e])
            for i in range(7):
                for j in range(7):
                    K[C[e][i]][C[e][j]] = K[C[e][i]][C[e][j]] + K_e[i][j]
                    
                    
        #-APPLY BC -->
        #--Reduced stiffness matrix
        K_red = np.copy(K)
        K_red = np.delete(K_red, BC, axis=0)
        K_red = np.delete(K_red, BC, axis=1)
        #--Reduced out of ballance force vector
        R_red = np.copy(R)
        R_red = np.delete(R_red, BC)
        #SOLVE REDUCED SYSTEM
        #-SOLVE FOR AN INCREMENT ddD IN dD, I.E. dD -> dD + ddD
        ddD_red = np.linalg.solve(K_red, R_red)
        
        #-UPDATE TOTAL DISPLACEMENT
        ddD = np.copy(ddD_red)
        for i in BC:
            ddD = np.insert(ddD, i, 0.)
        dD = dD + ddD
        D[inc+1] = D[inc]+dD
        
        
        #-UPDATE STRAIN (INCREMENT) AND THEN STRESS AT THE CORRESPONDING STRESS AND MATERIAL STIFFNESS 
        for e in range(Nel):
            crossSection_flag = elements[e][1]
            crossSection_el = crossSection_array[crossSection_flag]
            Nc_layers =  2*int(crossSections[crossSection_flag,2]) + 1
            crossSection_zc = crossSection_el['zc']
            Ns_layers =  int(crossSections[crossSection_flag,3])
            dD_e=np.zeros(7, dtype=float)
            for k in range(7):
                dD_e[k]=dD[C[e][k]]
            
            deps_allintPt = fu.strainCalc_concrete(crossSection_el, elCoord[e], dD_e)
            for i in range(2):
                for j in range(Nc_layers):
                    deps = deps_allintPt [i][j]
                    eps_c[e][inc+1][i][j] = eps_c[e][inc][i][j] + deps
                    predef = predef_array[inc][e][0] + crossSection_zc[j]*predef_array[inc][e][1]
                    dpredef = predef_array[inc+1][e][0] + crossSection_zc[j]*predef_array[inc+1][e][1]-predef
                    #-INTERNAL FORCES (STRESS UPDATE) 
                    #--CALL MATERIAL SUBROUTINES TO UPDATE STATE VARIABLES, STRESS AND MATERIAL STIFFNESS
                    sig_c[e][inc+1][i][j], St_c[e][i][j], stateVars_c[e][inc+1][i][j] = mat.matCalc_concrete(eps_c[e][inc][i][j], deps, time, dtime, predef, dpredef, matProps_c, sig_c[e][inc][i][j], St_c[e][i][j], stateVars_c[e][inc][i][j])
            
            deps_allintPt = fu.strainCalc_steel(crossSection_el, elCoord[e], dD_e)
            for i in range(2):
                for j in range(Ns_layers):
                    deps = deps_allintPt [i][j]
                    eps_s[e][inc+1][i][j] = eps_s[e][inc][i][j] + deps
        
                    #-INTERNAL FORCES (STRESS UPDATE) 
                    #--CALL MATERIAL SUBROUTINES TO UPDATE STATE VARIABLES, STRESS AND MATERIAL STIFFNESS
                    sig_s[e][inc+1][i][j], St_s[e][i][j], stateVars_s[e][inc+1][i][j] = mat.matCalc_steel(eps_s[e][inc][i][j], deps, time, dtime, predef, dpredef, matProps_s, sig_s[e][inc][i][j], St_s[e][i][j], stateVars_s[e][inc][i][j])

        #-UPDATE INTERNAL FORCE VECTOR
        Rint = np.zeros(Ndof, dtype=float)
        for e in range(Nel):
            crossSection_flag = elements[e][1]
            crossSection_el = crossSection_array[crossSection_flag]
            Rint_el = fu.elInternalForce(crossSection_el,elCoord[e], sig_c[e][inc+1], sig_s[e][inc+1])
            for i in range(7):
                Rint[C[e][i]] = Rint[C[e][i]] + Rint_el[i]

                
        #-OUT OF BALANCE FORCE, R  
        R = Rext - Rint
        #--Reduced out of ballance force vector
        R_red = np.copy(R)
        R_red = np.delete(R_red, BC)
        #-CONVERGENCE CHECK
        R_red_max = np.max( np.abs(R_red) )
        Rext_length2 = np.dot(Rext,Rext)
        if Rext_length2 > 0.:
            norm = np.sqrt(np.dot(R_red,R_red))/np.sqrt(Rext_length2)
        else:
            norm = R_red_max
        if R_red_max < 0.01:
        #if norm < 1E-7:
            string = 'CONVERGED INCREMENT ' + str(inc) + ' ON {:d} ITERATIONS'.format(itr)
            converegedIncrement = True
            lastConvergedIncrement = inc
            print(string)
            break
        else:
            print('MAXIMUM OUT OF BALLANCE FORCE IS: ', R_red_max)
    if converegedIncrement == False:
        print('ANALYSIS STOPPED!')
        print('NO CONVERGENCE ACHIEVED ON INCREMENT {:d}.'.format(inc))
        #break
        print('ANALYSIS WILL STILL CONTINUE.')
        
#-END SOLVING------------------------------------------------------------------    
    
    
#-PLOT AND SAVE SOME RESULTS IN RESULTS FOLDER
time_frame=lastConvergedIncrement+1

plt.close('all')
post.plot_disp(elCoord, C, D[Ninc_step1],D[time_frame], analysis_name)
post.plot_sectionForces_2(crossSection_array, elements, elCoord,sig_c,sig_s, analysis_name, time_frame)

#SAVE RESULTS
np.savez('RESULTS/{}'.format(analysis_name), elCoord=elCoord, sig_c=sig_c, stateVars_c=stateVars_c, sig_s=sig_s, stateVars_s=stateVars_s)