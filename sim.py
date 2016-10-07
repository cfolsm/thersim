#!/bin/bash/python


#-------------------------------------------------------------------------------------------------------------------------
#  Thermal Electric Simulation v0.1
#  Christopher Folmar
#  2016-09-06 
#  
#  All units are in SI standards, ie m, K, W ext...
#-------------------------------------------------------------------------------------------------------------------------


import sys
import cmath
import argparse
import csv
import math
import numpy as np
from pylab import plot,show
import scipy.fftpack
import matplotlib.pyplot as plt

sub = 'pvdf'
sub_thick = 1e-5
cond = 'ag'
sink = 'air'
debug = True 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CNT film characteristics 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha_var = np.linspace(10,150,1000)  
alpha = {'brd_champ': 110.5e-6, 'char_champ': 60.0e-6, 'var': alpha_var }  # Seebeck coefficient [V/K]
    # 191.3e-6 = 1500 uW/m*K^2, 184.8e-6 = 1400, 178e-6 = 1300, 171.1 = 1200
    # 163.8 = 1100, 156.2 = 1000, 148.2 = 900,139.7 = 800
    # 130.7 = 700, 121 = 600, 110.5 = 500, 98.8 = 400

film_sigma = 41000      # electrical conductivity of film [S/m]
film_sigma_var = np.linspace(1000,500000,1000)
film_thick = 1e-4  # film thickness [m]
film_thick_var = np.linspace(1e-9,1e-2,10000000)
T_h = 303            # Hot side temp [K]
T_c= 293           # Cold side temp [K]
sc = 'brd_champ'
#film_thick = film_thick_var
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Electrical resistivity constants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

film_rho = 1/film_sigma
rho = {'ag': 1.59e-8, 'cu': 1.68e-8, 'ag_ink': 33.8e-8, 'au': 2.2e-8, 'film': 1/41000 }
Ag_rho =  1.59e-8
Cu_rho = 1.68e-8
Ag_ink_rho = 33.8e-8
Au_rho = 2.2e-8
Ag_CNT_cr_var = np.linspace(1.59e-10,1e-2,10000000)
Au_CNT_cr_var = np.linspace(1.59e-8,1e-2,1000)
Ag_CNT_cr = 1e-6

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Thermal conductivity constants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

kappa = {'pvdf': 0.19, 'pet': 0.2, 'ptfe': 0.25, 'kapton': .00524, 'ag': 429, 'cu': 401, 'au': 310, 'ag_ink': 90, 'film': .1, 'air': .024} # thermal conductivities [W/m*K]
film_kappa = 1          # thermal conductivity of TEG (W/m*K)  
foam_kappa = 0.03       # Probably quite close to air which is 0.025 W/m*K
#Ag_kappa = 90          # Roughly 20% of bulk Ag at 429 W/m*K



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now we define the area of the entire structure, as well as the needed physical dimensions of the support
# structure.  We also define the physical parameters of our heat sinks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub_thick = 2.54e-5
sub_area = 1e-4

sink_thick = 1
sink_area = 1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We now wish to calculate the number of "modules" we can fit in a given area for both
# design layouts.  We will normalize these to 1cm^2.
# TF = Thin Film layout
# VS = Vertical stack layout
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#TF first
tf_area = 1.56e-7   #area of individual n and p module [m]
tf_area_var = np.linspace(1e-15,1e-5,10000000)
tf_area = tf_area_var
n_teg = sub_area/tf_area
n_teg_var = sub_area/tf_area_var

n_teg = n_teg_var

print('number of n-p modules',n_teg)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Electrical conductor parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cond = 'ag'
cond_thick = 1000e-9 #np.linspace(1e-9,1e-4,100000)  #thickness of conductor layer [m]
cond_width = tf_area**(1/2)  # [m]
cond_area = cond_thick*cond_width # [m^2]
cond_len = cond_width*2   # [m]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Thermal resistances
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Th_cond = (cond_thick)/(kappa[cond]*cond_area)*n_teg		#Thermal resistance of conductor [K/W]
print('conductor area ', cond_area)
print('conductor len ' , cond_len)
print('conductor thermal resistance = ', Th_cond)

Th_sub = (sub_thick)/(kappa[sub]*sub_area)
print('substrate thermal resistance = ', Th_sub)

Th_sink = (sink_thick)/(kappa[sink]*sink_area)

Th_Sl = (film_thick)/(kappa['film']*tf_area)*n_teg	#Thermal resistance of an individual n or p module [K/W]
print('film area ', tf_area)
print('film height ' , film_thick)

Th_Sl_var = 1/(kappa['film']*tf_area_var)

thermcont = np.linspace(1e-5,1,1000)			#Thermal contact resistance [m^2*K/W]
Th_cond = {'const': Th_cond , 'var': thermcont}	   	
Th_con = (Th_cond['const']+2*Th_sub) #2*Th_sink				#Thermal resistance of conductor layer [m^2*K/W]
Th_TEG = n_teg*(Th_con+Th_Sl) + 2*Th_sub




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TEG Electrical resistance 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

R_teg = (4*Ag_CNT_cr + (film_thick/tf_area)*film_rho + (cond_len/cond_area)*rho[cond])*n_teg  #Resistance of 1 TEG "module" [Ohm]

#R_ld = np.linspace(0.0001,2.9*R_teg,1000000)
#R_ld = 2.5*R_teg


alpha = n_teg*alpha[sc]
#Th_Sl = Th_Sl*n_teg
#Th_con = 0*Th_con*n_teg


if debug:
	print('resistance = ',R_teg)
#	print('Th_TEG = ', Th_TEG)
	print('Th_con = ', Th_con)
	print('Th_Sl = ' , Th_Sl)
	print('T_c = ', T_c)
	print('T_h = ', T_h)


print('seebeck = ', alpha)

#alpha = .0531876
#R_teg = 1.6
#Th_Sl = 1.498
#Th_con = np.linspace(.1,1.5,1000)
T_h = 770
T_c = 700
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  voltage and current calculation for a p-n structure 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


A = (Th_con + Th_Sl)*((alpha**2)*Th_con*Th_Sl*(T_h+T_c) + R_teg*(2*Th_con + Th_Sl))
#print('A = ',A)

B = ((Th_Sl + 2*Th_con)**2)*((alpha**2)*Th_con*Th_Sl*(T_h+T_c) + R_teg*(Th_Sl +2*Th_con))**2 - (alpha**4)*(Th_con**2)*(Th_Sl**3)*((T_h-T_c)**2)*(Th_Sl + 2*Th_con)
#print('B = ',B)

I_mpp = (A - np.sqrt(B))/((alpha**3)*(Th_con**2)*(Th_Sl**2)*(T_h-T_c))
#print('I_mpp = ',I_mpp)

P_max = (I_mpp/((alpha**2)*(Th_con**2)*(I_mpp**2) - Th_Sl - 2*Th_con))*(((alpha**2)*Th_con*Th_Sl*(T_h+T_c) + R_teg*(Th_Sl + 2*Th_con))*I_mpp - alpha*Th_Sl*(T_h-T_c))

#I_mpp = ((Th_Sl/(Th_Sl+2*Th_con))*alpha*(T_h-T_c))/(2*(R_teg+(((alpha**2)*Th_con*Th_Sl)/(Th_Sl+2*Th_con))*(T_h+T_c)))

#P_max = (((Th_Sl/(Th_Sl+2*Th_con))*alpha*(T_h-T_c))**2)/(4*(R_teg+(((alpha**2)*Th_con*Th_Sl)/(Th_Sl+2*Th_con))*(T_h+T_c)))

V_teg = P_max/I_mpp

#V_teg = (alpha[sc]*T_h*R_ld - alpha[sc]*T_c*R_ld)/(R_teg+R_ld) #no thermal contact resitance

#V_teg = (2*3**(1/3)*alpha**4*Th_con**3*(T_c+T_h)*R_teg*R_ld**2*Th_Sl**2 + 2*3**(1/3)*alpha**2*Th_con**2*R_teg*R_ld**2*(R_teg+R_ld)*Th_Sl*(2*Th_con+Th_Sl) - 2**(1/3)*((9*alpha**5*Th_con**4*(T_h - T_c)*R_teg**2*R_ld**3*Th_Sl**3 + 3**(1/2)*(alpha**6*Th_con**6*R_teg**3*R_ld**6*Th_Sl**3*(27*alpha**4*Th_con**2*(T_c-T_h)**2*R_teg*Th_Sl**3 - 4*(2*Th_con*(R_teg+R_ld) + alpha**2*Th_con*(T_c+T_h)*Th_Sl + (R_teg+R_ld)*Th_Sl)**3)))**(1/2))**(2/3))/((6**(2/3)*alpha**2*Th_con**2*R_teg*Th_Sl*(9*alpha**5*Th_con**4*(T_h-T_c)*R_teg**2*R_ld**3*Th_Sl**3 + 3**(1/2)*(alpha**6*Th_con**6*R_teg**3*R_ld**6*Th_Sl**3*(27*alpha**4*Th_con**2*(T_c-T_h)**2*R_teg*Th_Sl**3 - 4*(2*Th_con*(R_teg+R_ld) + alpha**2*Th_con*(T_c+T_h)*Th_Sl + (R_teg+R_ld)*Th_Sl)**3)))**(1/2))**(1/3))

#V_teg = -alpha*T_c*R_ld**3*Th_Sl + alpha*T_h*R_ld**3*Th_Sl + (-2*Th_con*R_teg*R_ld**2 - 2*Th_con*R_ld**3 - alpha**2*T_c*Th_con*R_ld**2*Th_Sl - alpha**2*Th_con*T_h*R_ld**2*Th_Sl - R_teg*R_ld**2*Th_Sl - R_ld**3*Th_Sl) + alpha**2*Th_con**2*R_teg*Th_Sl

#V_teg = *R_ld*((alpha[sc]*n_teg - 1)*T_c + (alpha[sc]*n_teg + 1)*T_h)**(1/3)/((alpha[sc]*n_teg)**(1/3)*(Th_con['const']*n_teg)**(2/3)*(R_teg*n_teg)**(1/3))/10

#V_teg = math.sqrt(pre_fact)*R_ld*((alpha[sc]*n_teg - 1)*T_c + (alpha[sc]*n_teg + 1)*T_h)**(1/3)/((alpha[sc]*n_teg)**(1/3)*(Th_con['const']*n_teg)**(2/3)*(R_teg*n_teg)**(1/3))/10

#neg = alpha[sc]**2*T_c*R_ld**2*Th_Sl*Th_con+alpha[sc]**2*T_c*R_ld**2*Th_Sl*Th_con - (R_teg*R_ld**2*Th_Sl + 2*R_teg*R_ld**2+Th_Sl*R_ld**3+2*R_ld**3*Th_con)

#V_teg = (1/(3*(2**(1/2)*alpha[sc]**2*Th_Sl*Th_con**2)))*((27*alpha[sc]**5*T_c*R_ld**2*Th_Sl**3*Th_con**4 - 27*alpha[sc]**5*T_h*R_ld**2*Th_Sl**3*Th_con**4 + ((27*alpha[sc]**5*T_c*R_ld**2*Th_Sl**3*Th_con**4 + 27*alpha[sc]**5*T_h*R_ld**2*Th_Sl**3*Th_con**4)**2 + 108*alpha[sc]**6*Th_Sl**3*Th_con**6*(-alpha[sc]**2*T_c*R_ld*Th_Sl*Th_con + alpha[sc]**2*T_h*R_ld*Th_Sl*Th_con + R_teg*R_ld*Th_Sl + 2*R_teg*R_ld*Th_con + R_ld**2*Th_con)**3)**(1/2)))**(1/3)- ( 2**(1/2)*(alpha[sc]**2*(-T_c)*R_ld*Th_Sl*Th_con - alpha[sc]**2*T_h*R_ld*Th_Sl*Th_con - R_teg*R_ld*Th_Sl - 2*R_teg*R_ld*Th_con - R_ld**2*Th_Sl - 2*R_ld**2*Th_con))/((27*alpha[sc]**5*T_c*R_ld**2*Th_Sl**3*Th_con**4 - 27*alpha[sc]**5*T_h*R_ld**2*Th_Sl**3*Th_con**4 + ((27*alpha[sc]**5*T_c*R_ld**2*Th_Sl**3*Th_con**4 + 27*alpha[sc]**5*T_h*R_ld**2*Th_Sl**3*Th_con**4)**2 + 108*alpha[sc]**6*Th_Sl**3*Th_con**6*(alpha[sc]**2*T_c*R_ld*Th_Sl*Th_con + alpha[sc]**2*T_h*R_ld*Th_Sl*Th_con + R_teg*R_ld*Th_Sl + 2*R_teg*R_ld*Th_con + R_ld**2*Th_Sl + 2*R_ld**2*Th_con)**3)**(1/2))**(1/3))

#print(V_teg)
#V_teg = V_teg.real
#V_teg = ((108*alpha[sc]**5*T_c*R_teg**2*R_ld**3*Th_Sl**3*Th_con**4 + 108*alpha[sc]**5*T_c*R_teg*R_ld**4*Th_Sl**3*Th_con**4+27*alpha[sc]**5*T_c*R_ld**5*Th_Sl**3*Th_con**4 - 108*alpha[sc]**5*T_h*R_teg**2*R_ld**3*Th_Sl**3*Th_con**4 - 108*alpha[sc]**5*T_h*R_teg*R_ld**4*Th_Sl**3*Th_con**4 - 27*alpha[sc]**5*T_h*R_ld**5*Th_Sl**3*Th_con**4 + ((108*alpha[sc]**5*T_c*R_teg**2*R_ld**3*Th_Sl**3*Th_con**4 + 108*alpha[sc]**5*T_c*R_teg*R_ld**4*Th_Sl**3*Th_con**4 + 27*alpha[sc]**5*T_c*R_ld**5*Th_Sl**3*Th_con**4 - 108*alpha[sc]**5*T_h*R_teg**2*R_ld**3*Th_Sl**3*Th_con**4 - 108*alpha[sc]**5*T_h*R_teg*R_ld**4*Th_Sl**3*Th_con**4 - 27*alpha[sc]**5*T_h*R_ld**5*Th_Sl**3*Th_con**4)**2 + 108*(2*alpha[sc]**2*R_teg*Th_Sl*Th_con**2 + alpha[sc]**2*R_ld*Th_Sl*Th_con**2)**3*(alpha[sc]**2*T_c*R_ld**2*Th_Sl*Th_con + alpha[sc]**2*T_h*R_ld**2*Th_Sl*Th_con - R_teg*R_ld**2*Th_Sl - 2*R_teg*R_ld**2*Th_con + R_ld**3*(-Th_Sl) - 2*R_ld**3*Th_con)**3)**(1/2))**(1/3))/(3*2**(1/3)*(2*alpha[sc]**2*R_teg*Th_Sl*Th_con**2 + alpha[sc]**2*R_ld*Th_Sl*Th_con**2)) - (2**(1/3)*(alpha[sc]**2*T_c*R_ld**2*Th_Sl*Th_con + alpha[sc]**2*T_h*R_ld**2*Th_Sl*Th_con - R_teg*R_ld**2*Th_Sl - 2*R_teg*R_ld**2*Th_con + R_ld**3*(-Th_Sl) - 2*R_ld**3*Th_con))/((108*alpha[sc]**5*T_c*R_teg**2*R_ld**3*Th_Sl**3*Th_con**4 + 108*alpha[sc]*T_c*R_teg*R_ld**4*Th_Sl**3*Th_con**4 + 27*alpha[sc]**5*T_c*R_ld**5*Th_Sl**3*Th_con**4 + ((108*alpha[sc]**5*T_c*R_teg**2*R_ld**3*Th_Sl**3*Th_con**4 + 108*alpha[sc]**5*T_c*R_teg*R_ld**4*Th_Sl**3*Th_con**4 + 27*alpha[sc]**5*T_c*R_ld**5*Th_Sl**3*Th_con**4 -108*alpha[sc]**5*T_h*R_teg**2*R_ld**3*Th_Sl**3*Th_con**4 - 108*alpha[sc]**5*T_h*R_teg*R_ld**4*Th_Sl**3*Th_con**4 - 27*alpha[sc]**5*T_h*R_ld**5*Th_Sl**3*Th_con**4)**2 + 108*(2*alpha[sc]**2*R_teg*Th_Sl*Th_con**2 + alpha[sc]**2*R_ld*Th_Sl*Th_con**2)**3*(alpha[sc]*T_c*R_ld**2*Th_Sl*Th_con + alpha[sc]**2*T_h*R_ld**2*Th_Sl*Th_con - R_teg*R_ld**2*Th_Sl - 2*R_teg*R_ld**2*Th_con + R_ld**3*(-Th_Sl) - 2*R_ld**3*Th_con)**3)**(1/2))**(1/3))

if debug:   
	print('V_mpp = ', V_teg)
	print('I_mpp = ', I_mpp)


print('power = ' ,P_max)
#pre_fact = (n_teg**2*alpha[sc]**2*Th_Sl*R_ld)/(R_ld**2+R_teg**2+R_ld*R_teg)

#numerator = (T_h-T_c - n_teg*alpha[sc]*(V_teg/R_ld)*Th_con*((V_teg**2/R_ld**2)*Th_con*R_teg + T_h + T_c))**2
#numerator = (T_h-T_c)**2 - alpha[sc]**2*n_teg**3*(V_teg**2/R_ld**2)*(Th_con)**2*((T_h-T_c)**2 + n_teg**2*(V_teg**4/R_ld**4)*(Th_con)**2*R_teg**2 + (T_h-T_c)**2*n_teg**2*(V_teg**2/R_ld**2)*Th_con*R_teg) - (T_h-T_c)**2*(n_teg**4*alpha[sc]*(V_teg/R_ld)*Th_con*((V_teg**2/R_ld**2)*Th_con*R_teg+T_h+T_c))

#divisor = ((2*Th_con + Th_Sl - (alpha[sc]*n_teg)**2*((V_teg**2/R_ld**2)*(Th_con**2*Th_Sl))))**2
#divisor = 2*Th_con**2*n_teg + Th_Sl**2*n_teg + 2*Th_con*Th_Sl*n_teg**2 + n_teg*2*alpha[sc]**2*(V_teg**2/R_ld**2)*(n_teg**8*alpha[sc]**2*(V_teg**2/R_ld**2)*Th_con**4*Th_Sl**2 - 2*Th_con**2*Th_Sl**2*n_teg**3 + Th_con**2*Th_Sl**2*n_teg**2)

#P1_teg = ((2*3**(1/3)*alpha**4*Th_con**3*(T_c+T_h)*R_teg*R_ld**2*Th_Sl**2 + 2*3**(1/3)*alpha**2*Th_con**2*R_teg*R_ld**2*(R_teg+R_ld)*Th_Sl*(2*Th_con+Th_Sl) - 2**(1/3)*(9*alpha**5*Th_con**4*(T_h-T_c)*R_teg**2*R_ld**3*Th_Sl**3 + 3**(1/2)*(alpha**6*Th_con**6*R_teg**3*R_ld**6*Th_Sl**3*(27*alpha**4*Th_con**2*(T_c-T_h)**2*R_teg*Th_Sl**3 - 4*(2*Th_con*(R_teg+R_ld) + alpha**2*Th_con*(T_c+T_h)*Th_Sl)**3)))**(1/2)**(2/3))**2/(6*6**(1/3)*alpha**4*Th_con**4*R_teg**2*R_ld*Th_Sl**2*(9*alpha**5*Th_con**4*(T_h-T_c)*R_teg**2*R_ld**3*Th_Sl**3 + 3**(1/2)*(alpha**6*Th_con**6*R_teg**3*R_ld**6*Th_Sl**3*(27*alpha**4*Th_con**2*(T_c-T_h)**2*R_teg*Th_Sl**3 - 4*(2*Th_con*(R_teg+R_ld) + alpha**2*Th_con*(T_c+T_h)*Th_Sl + (R_teg+R_ld)*Th_Sl)**3)))**(1/2)**(2/3))


#P = (V_teg**2)/R_ld
#print(P1_teg)

#P_teg = (1/R_ld)*pre_fact*(numerator/divisor)
#print('max R_ld = ', min(P))

'''
cond_height = 100e-9

v_m1 = alpha*n_teg*dT_TEAg/(1+2*(kappa['film']/kappa[cond])*(1e-9/film_height_var))
v_m2 = alpha*n_teg*dT_TEAg/(1+2*(kappa['film']/kappa[cond])*(5e-9/film_height_var))
v_m3 = alpha*n_teg*dT_TEAg/(1+2*(kappa['film']/kappa[cond])*(10e-9/film_height_var))
v_m4 = alpha*n_teg*dT_TEAg/(1+2*(kappa['film']/kappa[cond])*(cond_height/film_height_var))
v_m5 = alpha*n_teg*dT_TEAg/(1+2*(kappa['film']/kappa[cond])*(1000e-9/film_height_var))

i_m1 = alpha*(1/n_teg)*dT_TEAg/(2*rho['film']*((2*rho[cond]/rho['film'])+1)*(1+2*(kappa['film']/kappa[cond])*(1e-9/film_height_var)))
i_m2 = alpha*(1/n_teg)*dT_TEAg/(2*rho['film']*((2*rho[cond]/rho['film'])+1)*(1+2*(kappa['film']/kappa[cond])*(5e-9/film_height_var)))
i_m3 = alpha*(1/n_teg)*dT_TEAg/(2*rho['film']*((2*rho[cond]/rho['film'])+1)*(1+2*(kappa['film']/kappa[cond])*(10e-9/film_height_var)))
i_m4 = alpha*(1/n_teg)*dT_TEAg/(2*rho['film']*((2*rho[cond]/rho['film'])+1)*(1+2*(kappa['film']/kappa[cond])*(cond_height/film_height_var)))
i_m5 = alpha*(1/n_teg)*dT_TEAg/(2*rho['film']*((2*rho[cond]/rho['film'])+1)*(1+2*(kappa['film']/kappa[cond])*(1000e-9/film_height_var)))

#Max power achievable given heat exchange resistances
#P_max = pow(dT,2)*nu_reduced/(4*T_hot*(Rth_external)*Area_module) # Units should be W/cm^2
#Actual power given exact structure, assumes no electrical loses in Ag
#P_act = Q_TE*nu/Area_module # W/cm^2

p1 = v_m1*i_m1
p2 = v_m2*i_m2
p3 = v_m3*i_m3
p4 = v_m4*i_m4
p5 = v_m5*5

'''
fig, ax = plt.subplots()
#ax.plot(xf, 2.0/L * np.abs(A[:L/2]))
#ax.plot(film_height_var,p1)
#ax.plot(film_height_var,p2)
#ax.plot(film_height_var,p3)
#ax.plot(film_height_var,p4)
ax.plot(n_teg,P_max)
ax.set_xlabel('film thickness')
ax.set_ylabel('Power output [W]')
ax.set_xscale('log')

plt.show()


