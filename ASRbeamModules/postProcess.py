# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 14:56:08 2020

@author: simkon
"""
import numpy as np
import ASRbeamModules.functions as fu
import matplotlib.pyplot as plt
import matplotlib as mpl
import ASRbeamModules.my_plot as my_plot

#INPUT
latex_font_size = 10
latex_textwidth_pt = 345
latex_fraction_of_textwidth = 1.

plt.close('all')
nice_fonts = {
        # Use LaTeX to write all text
        #"pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
        "text.usetex": False,                # use LaTeX to write all text
        # "font.family": "serif",
        # "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
        # "font.sans-serif": [],
        # "font.monospace": [],
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": latex_font_size,
        "font.size": latex_font_size,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": latex_font_size,
        "xtick.labelsize": latex_font_size,
        "ytick.labelsize": latex_font_size,
        'lines.linewidth': 0.75,
        'lines.markeredgewidth': 0.75,
        'text.latex.preamble': r'\usepackage{txfonts}'
}
mpl.rcParams.update(nice_fonts)



#-PLOTTING

def plot_disp(elCoord, C, D1, D2, result_file_name):
    fig, ax  = plt.subplots() 
    Nel = np.shape(elCoord)[0]
    scale=1
    init_X_coord =np.array([])
    init_Z_coord =np.array([])
    disp1_X = np.array([])
    disp1_Z = np.array([])
    disp2_X = np.array([])
    disp2_Z = np.array([])
    for e in range(Nel):
        Du_e = np.array([D1[C[e][i]] for i in range(2)])
        Dw_e = np.array([D1[C[e][i]] for i in range(2,6)])
        
        Le = np.sqrt( (elCoord[e][2]-elCoord[e][0])**2. + (elCoord[e][3]-elCoord[e][1])**2.)
        s_theta = -(elCoord[e][3]-elCoord[e][1])/Le
        theta = np.arcsin(s_theta)
        
        d_e = np.matmul(fu.T(theta), np.append(Du_e,Dw_e))
        du_e = d_e[:2]
        dw_e = d_e[2:6]
        
        
        xe = np.linspace(0, Le,10)
        Xe = np.linspace(elCoord[e][0], elCoord[e][2], 10)
        Ze = np.linspace(elCoord[e][1], elCoord[e][3], 10)
        Ye = -Ze
        
        init_X_coord = np.append(init_X_coord,Xe)
        init_Z_coord = np.append(init_Z_coord,Ze)
    
        Ue = scale*np.array([fu.U(fu.u(i, Le, du_e), fu.w(i, Le, dw_e), theta) for i in xe])
        Ve = scale*np.array([fu.V(fu.u(i, Le, du_e), fu.w(i, Le, dw_e), theta) for i in xe])
        
        disp1_X = np.append(disp1_X,Ue)
        disp1_Z = np.append(disp1_Z,-Ve)

    
        ax.plot(np.add(Xe, Ue), np.add(Ye, Ve),  color='k', ls='-.' )
        ax.plot(Xe,Ye, 'k--')
        
    for e in range(Nel):
        Du_e = np.array([D2[C[e][i]] for i in range(2)])
        Dw_e = np.array([D2[C[e][i]] for i in range(2,6)])
        
        Le = np.sqrt( (elCoord[e][2]-elCoord[e][0])**2. + (elCoord[e][3]-elCoord[e][1])**2.)
        s_theta = -(elCoord[e][3]-elCoord[e][1])/Le
        theta = np.arcsin(s_theta)
        
        d_e = np.matmul(fu.T(theta), np.append(Du_e,Dw_e))
        du_e = d_e[:2]
        dw_e = d_e[2:6]
        
        
        xe = np.linspace(0, Le,10)
        Xe = np.linspace(elCoord[e][0], elCoord[e][2], 10)
        Ze = np.linspace(elCoord[e][1], elCoord[e][3], 10)
        Ye = -Ze
    
        Ue = scale*np.array([fu.U(fu.u(i, Le, du_e), fu.w(i, Le, dw_e), theta) for i in xe])
        Ve = scale*np.array([fu.V(fu.u(i, Le, du_e), fu.w(i, Le, dw_e), theta) for i in xe])
        
        disp2_X = np.append(disp2_X,Ue)
        disp2_Z = np.append(disp2_Z,-Ve)
    
        ax.plot(np.add(Xe, Ue), np.add(Ye, Ve), color='k', ls='-')

    
    ax.plot((0,0),(0,0), color='k', ls='-.', label=r'Dead load')
    ax.plot((0,0),(0,0), color='k', ls='-', label=r'Dead load + ASR expansion')
    ax.legend()
    ax.set_title(r'Deformed shape; displacements are scaled with {:.0f}'.format(scale))
    ax.set_frame_on(False)
    ax.axis('on')
    yticks = ax.get_yticks()
    ax.set_yticklabels(-yticks)
    ax.set_xlim(min(min(elCoord[:,0]),min(elCoord[:,2])), max(max(elCoord[:,0]),max(elCoord[:,2])))
    
    df = np.vstack((init_X_coord, init_Z_coord, disp1_X, disp1_Z, disp2_X, disp2_Z )).T
    np.savetxt("RESULTS/disp_{}.csv".format(result_file_name), df, delimiter=",")

    
    
def plot_sectionForces(crossSection_array, elements, elCoord, sig_c, sig_s, time_frame):
    Nel = np.shape(elCoord)[0]
    
    fig1=plt.figure()
    fig2=plt.figure()
    ax1=fig1.add_subplot(111)
    ax2=fig2.add_subplot(111)
    for e in range(Nel):
        crossSection_flag = elements[e][1]
        crossSection = crossSection_array[crossSection_flag]
        Le = np.sqrt( (elCoord[e][2]-elCoord[e][0])**2. + (elCoord[e][3]-elCoord[e][1])**2.)
        s_theta = -(elCoord[e][3]-elCoord[e][1])/Le
        theta = np.arcsin(s_theta)
        
        xe = np.linspace(0, Le,2)
        Xe = np.linspace(elCoord[e][0], elCoord[e][2], 2)
        Ze = np.linspace(elCoord[e][1], elCoord[e][3], 2)
        Ye = -Ze
        
        N_x1, M_x1 = fu.sectionForces(crossSection, sig_c[e][time_frame][0], sig_s[e][time_frame][0])
        N_x2, M_x2 = fu.sectionForces(crossSection, sig_c[e][time_frame][1], sig_s[e][time_frame][1])
        
        Me = 1/10**6*np.array([fu.epM(Le,  M_x1,  M_x2, i) for i in xe]) #kNm
        Ne = 1/10**3*np.array([fu.epM(Le,  N_x1,  N_x2, i) for i in xe]) #kN
        
        #-Transform the moment diagram for plotting purposes
        Me_X = np.sin(theta) * Me
        Me_Z = np.cos(theta) * Me
        Me_Y = -Me_Z
        
        Ne_X = np.sin(theta) * Ne
        Ne_Z = np.cos(theta) * Ne
        Ne_Y = -Ne_Z
        
        if e <62:
            color_of_plot = 'r'
        else:
            color_of_plot = 'b'
        scaler = 1
        ax1.plot(np.add(Xe, Me_X), np.add(Ye, Me_Y), color=color_of_plot)
        ax1.text(np.add(Xe, Me_X*scaler)[0], np.add(Ye, Me_Y*scaler)[0], r'$ %.1f$' % Me[0])
        ax1.plot(Xe,Ye, 'k--')
        ax2.plot(np.add(Xe, Ne_X), np.add(Ye, Ne_Y), color=color_of_plot)
        ax2.text(np.add(Xe, Ne_X*scaler)[0], np.add(Ye, Ne_Y*scaler)[0], r'$ %.1f$' % Ne[0])
        ax2.plot(Xe,Ye, 'k--')
        
def plot_sectionForces_2(crossSection_array, elements, elCoord, sig_c, sig_s, result_file_name, time_frame):
    Nel = np.shape(elCoord)[0]
    
    fig1=plt.figure()
    fig2=plt.figure()
    fig3=plt.figure()
    ax1=fig1.add_subplot(111)
    ax2=fig2.add_subplot(111)
    ax3=fig3.add_subplot(111)
    ax1.set_title(r'Bending moment')
    ax2.set_title(r'Axial force')
    ax3.set_title(r'Shear force')
    N_array = np.array([])
    M_array = np.array([])
    V_array = np.array([])
    Node_array =np.array([])
    for e in range(Nel):
        crossSection_flag = elements[e][1]
        crossSection = crossSection_array[crossSection_flag]
        Le = np.sqrt( (elCoord[e][2]-elCoord[e][0])**2. + (elCoord[e][3]-elCoord[e][1])**2.)
        s_theta = -(elCoord[e][3]-elCoord[e][1])/Le
        theta = np.arcsin(s_theta)
        Sy = 0.
        for j in range(len(crossSection['Ac'])):
            Sy = Sy + crossSection['zc'][2*j+1]*crossSection['Ac'][j]
        for j in range(len(crossSection['zs'])):
            Sy = Sy + crossSection['zs'][j]*crossSection['As'][j]
        
        #xe = np.linspace(0, Le,2)
        Xe = np.linspace(elCoord[e][0], elCoord[e][2], 2)
        Ze = np.linspace(elCoord[e][1], elCoord[e][3], 2)
        Ye = -Ze
        
        Rint = fu.elInternalForce(crossSection,elCoord[e], sig_c[e][time_frame], sig_s[e][time_frame])
        
        
        
        Ne = np.array([-Rint[0], Rint[1]]) 
        
        
        Me = 1/10**6*(np.array([Rint[3], -Rint[5]]) )#- Sy/(sum(crossSection['Ac'])+sum(crossSection['As'])) * Ne) #kNm
        Ne = 1/10**3*Ne #kN
        Ve = 1/10**3*np.array([-Rint[2], Rint[4]])
        
        N_array = np.append(N_array, Ne)
        M_array = np.append(M_array, Me)
        V_array = np.append(V_array, Ve)
        Node_array = np.append(Node_array, Xe)
        #-Transform the moment diagram for plotting purposes
        Me_X = np.sin(theta) * Me
        Me_Z = np.cos(theta) * Me
        Me_Y = -Me_Z
        
        Ne_X = np.sin(theta) * Ne
        Ne_Z = np.cos(theta) * Ne
        Ne_Y = -Ne_Z
        
        Ve_X = np.sin(theta) * Ve
        Ve_Z = np.cos(theta) * Ve
        Ve_Y = -Ve_Z
        
        scaler = 50
        ax1.plot(np.add(Xe, Me_X*scaler), np.add(Ye, Me_Y*scaler), color='0.85')
        ax1.set_xlabel(r'$x$ [mm]')
        ax1.set_ylabel(r'$M$ [kNm]')
        ax1.yaxis.set_ticks([])
        ax1.plot(Xe,Ye, 'k--')
        ax1.text(np.add(Xe, Me_X*scaler)[0], np.add(Ye, Me_Y*scaler)[0], r'$ %.1f$' % Me[0])
        ax2.plot(np.add(Xe, Ne_X), np.add(Ye, Ne_Y), color='0.85')
        ax2.plot(Xe,Ye, 'k--')
        ax3.plot(np.add(Xe, Ve_X), np.add(Ye, Ve_Y), color='0.85')
        ax3.plot(Xe,Ye, 'k--')
    
    df = np.vstack((Node_array, N_array, M_array, V_array )).T

    
    np.savetxt("RESULTS/nodeForce_{}.csv".format(result_file_name), df, delimiter=",")
    

def plot_openCracksAndSteelPlasticity(elCoord, sig_c, stateVars_c, stateVars_s, time_frame, result_file_name):
    Nel = int(np.shape(elCoord)[0]/2)
    NintPts = np.shape(stateVars_c[0][0])[1]
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex='col', figsize=my_plot.set_size(latex_textwidth_pt, latex_fraction_of_textwidth))
    axes_titles = (r'a) Inner beam',r'b) Outer beam')
    for k in range(2):
        ax = axes[k]
        ax.set_frame_on(False)
        ax.set_title(axes_titles[k], loc='left', fontsize=latex_font_size)
        
    
        for e in range(k*Nel,Nel+k*Nel):
            L = elCoord[e,2]-elCoord[e,0]
            x = np.array([elCoord[e,0] + L/2.-3.**0.5/6*L, elCoord[e,0] + L/2.+3.**0.5/6    *L])
            y = np.linspace(-1740/2, 1740/2,NintPts)
            for i in range(2):
                for j in range(NintPts):
                    if stateVars_c[e][time_frame][i][j][0] > 0. and sig_c[e][time_frame][i][j]>0.:
                        ax.scatter(x[i],y[j],  marker = (2, 2, 0.), color='k', s=10)
    
                if stateVars_s[e][time_frame][i][0][0]> 0.:
                    ax.scatter(x[i],-811,  marker = (2, 2, 90.), color='r', s=10)
                if stateVars_s[e][time_frame][i][1][0]> 0.:
                    ax.scatter(x[i], 811,  marker = (2, 2, 90.), color='r', s=10)
        ax.plot((0,0), (-1740/2, 1740/2), color='grey')       
        ax.plot((0,66250), (-1740/2, -1740/2), color='grey')
        ax.plot((0,66250), (1740/2, 1740/2), color='grey')   
        ax.plot((66250,66250), (-1740/2, 1740/2), color='grey')     
        ax.set_ylim(1500.,-1500)
        ax.set_xlim(-500.,66250+500)
        ax.scatter((0,22500,45000,66250), np.array([    1740/2,1740/2,1740/2,1740/2])+100, marker = '^', color = 'grey', label = r'Support')
        
        ax.set_ylabel(r'$z$ [mm]')
        ax.get_yaxis().set_visible(False)
    axes[0].get_xaxis().set_visible(False)
    ax.set_xlabel(r'$x$ [mm]')
    axes[0].scatter(-10000,0,  marker = (2, 2, 0.), color='k', label=r'Crack')
    axes[0].scatter(-10000,0,  marker = (2, 2, 90.), color='r', label=r'Yielding reinforcement')
    axes[0].legend(bbox_to_anchor=(0., (1500+1800/2)/3000, 1., (1500+1800/2)/3000), loc='lower left', ncol=3, mode="expand", borderaxespad=0., frameon=False)
    fig.tight_layout()
    fig.savefig('RESULTS/crackPattern_{}.pdf'.format(result_file_name))
    fig.savefig('RESULTS/crackPattern_{}.png'.format(result_file_name), dpi=500)

def plot_openCracksAndSteelPlasticity_innerBeam(elCoord, sig_c, stateVars_c, stateVars_s, time_frame, result_file_name):
    Nel = int(np.shape(elCoord)[0]/2)
    NintPts = np.shape(stateVars_c[0][0])[1]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.77376504773765, 2))


    ax.set_frame_on(False)

    for e in range(Nel):
        L = elCoord[e,2]-elCoord[e,0]
        x = np.array([elCoord[e,0] + L/2.-3.**0.5/6*L, elCoord[e,0] + L/2.+3.**0.5/6    *L])
        y = np.linspace(-1740/2, 1740/2,NintPts)
        for i in range(2):
            for j in range(NintPts):
                if stateVars_c[e][time_frame][i][j][0] > 0. and sig_c[e][time_frame][i][j]>0.:
                    ax.scatter(x[i],y[j],  marker = (2, 2, 0.), color='k', s=10)

            if stateVars_s[e][time_frame][i][0][0]> 0.:
                ax.scatter(x[i],-811,  marker = (2, 2, 90.), color='r', s=10)
            if stateVars_s[e][time_frame][i][1][0]> 0.:
                ax.scatter(x[i], 811,  marker = (2, 2, 90.), color='r', s=10)
    ax.plot((0,0), (-1740/2, 1740/2), color='grey')       
    ax.plot((0,66250), (-1740/2, -1740/2), color='grey')
    ax.plot((0,66250), (1740/2, 1740/2), color='grey')   
    ax.plot((66250,66250), (-1740/2, 1740/2), color='grey')     
    ax.set_ylim(1500.,-1500)
    ax.set_xlim(-500.,66250+500)
    ax.scatter((0,22500,45000,66250), np.array([    1740/2,1740/2,1740/2,1740/2])+100, marker = '^', color = 'grey', label = r'Support')
    
    ax.set_ylabel(r'$z$ [mm]')
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.set_xlabel(r'$x$ [mm]')
    ax.scatter(-10000,0,  marker = (2, 2, 0.), color='k', label=r'Crack')
    ax.scatter(-10000,0,  marker = (2, 2, 90.), color='r', label=r'Yielding reinforcement')
    ax.legend(bbox_to_anchor=(0., (1500+1800/2)/3000, 1., (1500+1800/2)/3000), loc='lower left', ncol=3, mode="expand", borderaxespad=0., frameon=False)
    ax.set_aspect(aspect=5)
    fig.tight_layout()
    fig.savefig('RESULTS/crackPattern_{}.pdf'.format(result_file_name), format='pdf', bbox_inches='tight')
    fig.savefig('RESULTS/crackPattern_{}.png'.format(result_file_name), dpi=500)



def plot_openCracksAndSteelPlasticity_innerBeam_2(elCoord, sig_c, stateVars_c, stateVars_s, time_frames, result_file_name):
    Nel = int(np.shape(elCoord)[0]/2)
    NintPts = np.shape(stateVars_c[0][0])[1]
    
    fig, ax = plt.subplots(nrows=len(time_frames), ncols=1, figsize=(4.77376504773765, 2))
    axes_titles = (r'a) After dead load', r'b) After ASR expansion')
    for q in range(len(time_frames)):
        time_frame = time_frames[q]
        ax[q].set_frame_on(False)
        ax[q].set_title(axes_titles[q], loc='left', fontsize=latex_font_size)
        for e in range(Nel):
            L = elCoord[e,2]-elCoord[e,0]
            x = np.array([elCoord[e,0] + L/2.-3.**0.5/6*L, elCoord[e,0] + L/2.+3.**0.5/6    *L])
            y = np.linspace(-1740/2, 1740/2,NintPts)
            for i in range(2):
                for j in range(NintPts):
                    if stateVars_c[e][time_frame][i][j][0] > 0. and sig_c[e][time_frame][i][j]>0.:
                        ax[q].scatter(x[i],y[j],  marker = (2, 2, 0.), color='k', s=10)
    
                if stateVars_s[e][time_frame][i][0][0]> 0.:
                    ax[q].scatter(x[i],-811,  marker = (2, 2, 90.), color='r', s=10)
                if stateVars_s[e][time_frame][i][1][0]> 0.:
                    ax[q].scatter(x[i], 811,  marker = (2, 2, 90.), color='r', s=10)
        ax[q].plot((0,0), (-1740/2, 1740/2), color='grey')       
        ax[q].plot((0,66250), (-1740/2, -1740/2), color='grey')
        ax[q].plot((0,66250), (1740/2, 1740/2), color='grey')   
        ax[q].plot((66250,66250), (-1740/2, 1740/2), color='grey')     
        ax[q].set_ylim(1500.,-1500)
        ax[q].set_xlim(-500.,66250+500)
        ax[q].scatter((0,22500,45000,66250), np.array([    1740/2,1740/2,1740/2,1740/2])+100, marker = '^', color = 'grey', label = r'Support')
        
        ax[q].set_ylabel(r'$z$ [mm]')
        ax[q].get_yaxis().set_visible(False)
        ax[q].get_xaxis().set_visible(False)
        ax[q].set_xlabel(r'$x$ [mm]')
        ax[q].set_aspect(aspect=5)
    
    ax[0].scatter(-10000,0,  marker = (2, 2, 0.), color='k', label=r'Crack')
    ax[0].scatter(-10000,0,  marker = (2, 2, 90.), color='r', label=r'Yielding reinforcement')
    ax[0].legend(bbox_to_anchor=(0., (1500+1800/2)/3000, 1., (1500+1800/2)/3000), loc='lower left', ncol=3, mode="expand", borderaxespad=0., frameon=False)
    fig.tight_layout()
    fig.savefig('RESULTS/crackPattern_{}.pdf'.format(result_file_name), format='pdf', bbox_inches='tight')
    fig.savefig('RESULTS/crackPattern_{}.png'.format(result_file_name), dpi=500)

#def crack_update(time_frame):
#    ax.clear()
#    ax.set_title(r'Time step {}'.format(time_frame))
#    ax.set_frame_on(False)
#    ax.plot((0,0), (-1740/2, 1740/2), color='grey')       
#    ax.plot((0,66250), (-1740/2, -1740/2), color='grey')
#    ax.plot((0,66250), (1740/2, 1740/2), color='grey')   
#    ax.plot((66250,66250), (-1740/2, 1740/2), color='grey')     
#    ax.set_ylim(4000.,-4000)
#    ax.set_xlim(-500.,66250+500)
#    ax.scatter((0,22500,45000,66250), np.array([    1740/2,1740/2,1740/2,1740/2])+100, marker = '^', color = 'grey', label = r'Support')
#    
#    ax.set_ylabel(r'$z$ [mm]')
#    ax.get_yaxis().set_visible(False)
#    
#    ax.get_xaxis().set_visible(False)
#    ax.set_xlabel(r'$x$ [mm]')
#    Nel = int(np.shape(elCoord)[0]/2)
#    NintPts = np.shape(stateVars_c[0][0])[1]
#    crack_list = np.array([[-10000,0]])
#    yield_list = np.array([[-10000,0]])
#
#    for e in range(0,Nel):
#        L = elCoord[e,2]-elCoord[e,0]
#        x = np.array([elCoord[e,0] + L/2.-3.**0.5/6*L, elCoord[e,0] + L/2.+3.**0.5/6    *L])
#        y = np.linspace(-1740/2, 1740/2,NintPts)
#        for i in range(2):
#            for j in range(NintPts):
#                if stateVars_c[e][time_frame][i][j][0] > 0. and sig_c[e][time_frame][i][j]>0.:
#                    crack_list = np.append(crack_list, np.array([[ x[i], y[j] ]]), axis=0)
#            
#            if stateVars_s[e][time_frame][i][0][0]> 0.:
#                yield_list = np.append(yield_list, np.array([[x[i], -811.]]), axis=0) 
#            if stateVars_s[e][time_frame][i][1][0]> 0.:
#                yield_list = np.append(yield_list, np.array([[x[i], 811.]]), axis=0) 
#            
#    ax.scatter(yield_list[:,0], yield_list[:,1],  marker = (2, 2, 90.), color='r')
#    ax.scatter(crack_list[:,0], crack_list[:,1],  marker = (2, 2, 0.), color='k')
#   
#
#
#crack_list = np.array([[-10000,0]])
#yield_list = np.array([[-10000,0]])
#fig, ax = plt.subplots()
#
#ax.set_frame_on(False)
#ax.plot((0,0), (-1740/2, 1740/2), color='grey')       
#ax.plot((0,66250), (-1740/2, -1740/2), color='grey')
#ax.plot((0,66250), (1740/2, 1740/2), color='grey')   
#ax.plot((66250,66250), (-1740/2, 1740/2), color='grey')     
#ax.set_ylim(4000.,-4000)
#ax.set_xlim(-500.,66250+500)
#ax.scatter((0,22500,45000,66250), np.array([    1740/2,1740/2,1740/2,1740/2])+100, marker = '^', color = 'grey', label = r'Support')
#
#ax.set_ylabel(r'$z$ [mm]')
#ax.get_yaxis().set_visible(False)
#
#ax.get_xaxis().set_visible(False)
#ax.set_xlabel(r'$x$ [mm]')
#
#ax.scatter(crack_list[:,0], crack_list[:,1],  marker = (2, 2, 0.), color='k')
#ax.scatter(yield_list[:,0], yield_list[:,1],  marker = (2, 2, 90.), color='r')
#
#animation = mpl.animation.FuncAnimation(fig, crack_update, frames=35, interval=2000, repeat_delay=3000)
#animation.save('cracks.avi', extra_args=['-vcodec', 'libx264'])
