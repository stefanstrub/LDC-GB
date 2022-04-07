from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint
import sys
from itertools import combinations
import itertools
from scipy import linalg
import mpmath
from subprocess import call
from mpl_toolkits.mplot3d import Axes3D
import subprocess
from scipy.optimize import fsolve
from scipy import interpolate
import matplotlib as mpl
from matplotlib import colors as ccolor
from matplotlib import cm
import re    
import csv
from scipy.integrate import solve_ivp
from matplotlib import rcParams


#------------------------------------------------------------------------------------------
#Units and conversions:
#------------------------------------------------------------------------------------------
#code units: Rsun, Msun, G=1, ...
c_SI        = 299792458.0       #m/s
M_sun_SI    = 1.989*(10.**30.)  #kg
R_sun_SI    = 695800000.        #m
AU_SI       = 149597871000.     #m 
G_new_SI    = 6.67*(10.**(-11.))
AU_U        = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U     = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
sec_year    = 31536000.
m_parsec    = 3.086*(10**16.)   #m
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
#CALC ydot:
#------------------------------------------------------------------------------------------
def func_NBODY_Ydot(Yin_arr, t):    #incl masses later on
    #reshape:
    Yin_nrobj_posvel    = np.reshape(Yin_arr,(nr_obj,6))
    #define:
    Ydot_nrobj_posvel   = np.zeros((nr_obj,6), dtype='d')
    PN_gamma    = 9.86719698246e-09
    for i in range(0,nr_obj):
        a_ij_tot    = np.array([0,0,0]) 
        pos_i       = Yin_nrobj_posvel[i,0:3]
        vel_i       = Yin_nrobj_posvel[i,3:6]
        for j in range(0,nr_obj):
            if (i != j):
                mi      = m0_arr[i]
                mj      = m0_arr[j]
                pos_j   = Yin_nrobj_posvel[j,0:3]
                vel_j   = Yin_nrobj_posvel[j,3:6]
                pos_ij  = pos_j - pos_i
                vel_ij  = vel_j - vel_i
                #print np.sqrt(np.sum(pos_ij[:]**2.)), i,j
                r_ij    = np.sqrt(np.sum(pos_ij**2.))
                v_ij    = np.sqrt(np.sum(vel_ij**2.))
                #Newtonian acc:
                a_ij        = (mj/(r_ij**2.))*(pos_ij/r_ij)
                #2.5PN acc:
                n_ij		= pos_ij/r_ij
                n_dot_vi	= np.dot(n_ij,vel_i)
                n_dot_vj	= np.dot(n_ij,vel_j)
                a25_ij      = -(PN_gamma**(5./2.))*(4./5.)*(mi*mj/(r_ij**3))*(vel_ij*(-(v_ij**2.) + 2.*(mi/r_ij) - 8.*(mj/r_ij)) + n_ij*(n_dot_vi-n_dot_vj)*(3.*(v_ij**2.) - 6.*(mi/r_ij) + (52./3.)*(mj/r_ij)))                                                        
                #total acc:                                                                    
                a_ij_tot    = a_ij_tot + a_ij + a25_ij
        a_i_tot = a_ij_tot
        Ydot_nrobj_posvel[i,0:3]    = np.ravel(vel_i)
        Ydot_nrobj_posvel[i,3:6]    = np.ravel(a_i_tot)
    Yout_arr    = np.ravel(Ydot_nrobj_posvel)
    return Yout_arr
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#calc dt:
#------------------------------------------------------------------------------------------
def func_dt(Yin_arr):    #incl masses later on
    #reshape:
    Yin_nrobj_posvel    = np.reshape(Yin_arr,(nr_obj,6))
    dt_fac              = 1e10
    for i in range(0,nr_obj):
        for j in range(0,nr_obj):
            if (i != j):
                pos_ij      = Yin_nrobj_posvel[i,0:3] - Yin_nrobj_posvel[j,0:3]
                vel_ij      = Yin_nrobj_posvel[i,3:6] - Yin_nrobj_posvel[j,3:6]
                r_ij        = np.sqrt(np.sum(pos_ij**2.))
                v_ij        = np.sqrt(np.sum(vel_ij**2.))
                a_ij        = m0_arr[j]/r_ij**2
                dt_rvra     = [r_ij/v_ij, np.sqrt(r_ij/a_ij)]
                min_dt_rvra = min(dt_rvra)
                if (min_dt_rvra <= dt_fac):
                    dt_fac = min_dt_rvra
    return dt_fac
#------------------------------------------------------------------------------------------


#-----------------------------------------------------------------
#input:
#-----------------------------------------------------------------
nr_obj  = 3
m0_arr  = np.zeros(nr_obj,      dtype='d') 
Y0_arr  = np.zeros(nr_obj*6,    dtype='d')
#-----------------------------------------------------------------
#make ICs:
#-----------------------------------------------------------------
m1      = 2.0   #in Msun
m2      = 1.0   #in Msun
m3      = 1   #in Msun
SMA_bin = 1   #in AU
#-----------------------------------------------------------------
#define:
#-----------------------------------------------------------------
m_bin    = m1 + m2
#-----------------------------------------------------------------
#Set b1 and b2 up in a CIRCULAR binary with COM in (0,0,0):
#-----------------------------------------------------------------
v_redmass           = np.sqrt(m_bin/SMA_bin)
b1_posxyz_binCM     = np.array([ (m2/m_bin)*SMA_bin,0,0])
b2_posxyz_binCM     = np.array([-(m1/m_bin)*SMA_bin,0,0])
b1_velxyz_binCM     = np.array([0, 1.0*(m2/m_bin)*v_redmass,0])
b2_velxyz_binCM     = np.array([0,-1.0*(m1/m_bin)*v_redmass,0])
#-----------------------------------------------------------------
#IC for b3:
#-----------------------------------------------------------------
b3_posxyz_binCM     = np.array([2.*SMA_bin,0,0])
b3_velxyz_binCM     = np.array([-0.1*v_redmass,0,0])
#-----------------------------------------------------------------
m0_arr[:] = np.array([m1,m2,m3])
Y0_arr[0*6:0*6 + 6] = np.append(b1_posxyz_binCM, b1_velxyz_binCM)
Y0_arr[1*6:1*6 + 6] = np.append(b2_posxyz_binCM, b2_velxyz_binCM)
Y0_arr[2*6:2*6 + 6] = np.append(b3_posxyz_binCM, b3_velxyz_binCM) 
#-----------------------------------------------------------------


#-----------------------------------------------------------------
#EVOLVE SYSTEM:
#-----------------------------------------------------------------
Nsteps      = 1000
T_orb       = 2.*np.pi*np.sqrt((SMA_bin**3.)/m_bin)
T_max       = 50.*T_orb
dt_scale    = 0.5

fig = plt.figure(figsize=(4, 4))
ax1  = fig.add_subplot(111)

pos_xyz_allt   = np.zeros((Nsteps,nr_obj,3), dtype='d')

Y_evolv     = Y0_arr
for sc in range(0,Nsteps):
    
    #variable dt:
    dt_step = func_dt(Y_evolv)
    #choose dt:
    dt_evolve   = dt_scale*dt_step 
    #dt_evolve   = dt_scale*T_orb
    
    #evolve system:
    Ydt_evolve  = odeint(func_NBODY_Ydot, Y_evolv, np.array([0.0, dt_evolve]), full_output=0, atol=1e-6, rtol=1e-6)[1]
    Y_evolv     = Ydt_evolve  
    
    #save for plotting: 
    pos_xyz_allt[sc,0,:]    = Y_evolv[0:3] 
    pos_xyz_allt[sc,1,:]    = Y_evolv[6:9] 
    pos_xyz_allt[sc,2,:]    = Y_evolv[12:15] 

#PLOT:
ax1.plot(pos_xyz_allt[:,0,0], pos_xyz_allt[:,0,1], marker = 'o', markersize = 0.1, color = 'black')
ax1.plot(pos_xyz_allt[:,1,0], pos_xyz_allt[:,1,1], marker = 'o', markersize = 0.1, color = 'red')
ax1.plot(pos_xyz_allt[:,2,0], pos_xyz_allt[:,2,1], marker = 'o', markersize = 0.1, color = 'blue')

ax1.set_xlim(-8*SMA_bin, 8*SMA_bin)
ax1.set_ylim(-8*SMA_bin, 8*SMA_bin)

plt.show()




















