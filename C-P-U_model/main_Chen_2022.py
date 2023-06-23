#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A model for the C-P-U evolution 310-285 Ma
Author: Shihan Li
"""
import numpy as np
from scipy.integrate import solve_ivp  
import timeit


import load_forcings as ffunc         # COPSE_reloaded forcings
import functions as fts               # functions for ODE

""" --------------------------------------------------------    Constant    ------------------------------------------------ """
k_CtoK = 273.15    # convert degrees C to Kelvin
k_c = 4.3280                 # new determines climate sensitivity to CO2
k_l = 7.4000                 # new determines temperature sensitivity to luminosity
k_planenhance = 0.15         # plant enhancement of weathering
f_oxw_a = 0.5                # oxidative weathering dependency on O2 concentration
k_logistic = 12              # new determines slope of logistic anoxia function, COPSE_reloaded, Clarkson et al. 2018
k_uptake = 0.5000            # new determines efficiency of nutrient uptake in anoxia function; COPSE_reloaded; Clarkson et al. 2018
# newp0 = 225.956              # new production (umol/kg)
newp0 = 117 * 2.2         
k_oxfrac = 0.9975            # updated initial oxic fraction
f_mocb_b = 2                 # marine organic carbon burial dependency on new production
k_oceanmass = 1.397e21       # ocean mass (kg)
""" --------------------------------------------------------  Modern Fluxes ------------------------------------------------ """

# Carbon, processes considered here: volcanic degassing, weathering (organic, carbonate, silicate), burial (continental organic, marine organic, carbonate)
k_ccdeg = 8e12      # volcanic degassing
k_silw = 8e12       # silicate weathering
k_carbw = 8e12      # carbonate weathering
k_mccb = k_silw + k_carbw # carbonate burial

k_oxidw = 5e12      # oxidative weatheirng + degassing 
k_locb = 2.5e12     # continental organic C burial
k_mocb = 2.5e12     # marine organic C burial

pco2atm0 = 280e-6      # pco2 corresponding to a0 (atm)
a0 = 3.193e18       # atmosphere-ocean co2
o0 = 3.7e19         # atmosphere-ocean o2

# P
p0 = 3.1e15         # ocean (phosphate) phosphorus
# n0 = 4.35e16        # ocean (reactive) nitrogen
k_fepb = 1e10       # fe-p burial
k_capb = 2*k_fepb      # ca-p burial
k_mopb = k_fepb      # organic-P burial
k_phosw = 4*k_fepb      # phosphate weathering 
   

# U
u0 = 1.85e13         # modern U in the ocean

u_riv0 = 4.79e7       # river input
u_hydro0 = 5.7e6      # hydrothermal output
u_anox0 = 6.2e6       # anoxic sink
u_sed0 = 3.6e7         
    
# isotope fractionations
d238u_riv = -0.29
u_frac_anox = 0.5
u_frac_sed = 0.0156
u_frac_hydro = 0.2










""" --------------------------------------------------------  Initial values (310 Ma) ------------------------------------------------ """
t_i = -310e6
t_Ma_i = t_i/1e6

pco2_i = 500e-6                      # after proxy data
pco2pal_i = pco2_i/pco2atm0     
a_i = a0 * np.sqrt(pco2pal_i)
temp_i = k_CtoK + 15 + k_c * np.log(pco2pal_i) + k_l * t_i/570e6   # luminosity used imiplicitly
                                                                   # corrected 0 degrees by the sea surface temperature record


o_i = 4.4e19         # after COPSE result, equivalent to 25% O2atm, o2_pal 1.19

po2pal_i = o_i/o0    # relative o2 concentration 


# after proxy
delta_ocn_i = 4.42
delta_u_i = -0.14

# initial degassing 


#### initial flux 


# forcing parameter, degassing, uplift driving erosion, effects of vegetation and lithology on weatherability
PG_i = ffunc.fPG(t_Ma_i)                        # Paleogeography effect on runoff/weathering, (Royer et al., 2014)
DEGASS_i = ffunc.fhaq_D(t_Ma_i)                 # metamorphic and volcanic degassing, inversion of sea-level curve. (Mills et al., 2017)
# DEGASS_i = 1.14                                # GEOCARBSURF 
# DEGASS_i = 1.  
UPLIFT_i = ffunc.fUPLIFT(t_Ma_i)                # tectonic uplift from sediment accumulation rates
# UPLIFT_i = 0.7743                               # Royer et al. 2014
EVO_i = ffunc.fEVO(t_Ma_i)                      # Land plant evolution
W_i = ffunc.fW(t_Ma_i)                          # Land plant enhancement of weathering
# W_i = 0.875                                   # GEOCARBSURF

# degassing
ccdeg_i = k_ccdeg * DEGASS_i

# silicate weathering, assumed equal to degassing for a steady state
silw_i = ccdeg_i

# biomass

VEG_i = fts.land_biota(temp_i-k_CtoK, pco2_i, po2pal_i, EVO_i) 

# weathering rates, activation energy as 62 kJ/mol
f_T_silw_i = np.exp(0.09 * (temp_i-288.15))
f_T_ap_i = np.exp(0.0507 * (temp_i - 288.15))  # 35 kJ/mol

# runoff temperature dependence
f_T_runoff_i = ((1 + 0.038 * (temp_i - 288.15))**0.65)


f_co2_abiotic_i = pco2pal_i ** 0.5 

f_co2_biotic = (2*pco2pal_i / (1+pco2pal_i))**0.4

# # o2 dependence?
VWmin_i = min(VEG_i * W_i, 1)
# f_silw_i = f_T_silw_i * f_T_runoff_i * ((1-VWmin_i) * k_planenhance * f_co2_abiotic_i + VEG_i * W_i)

# silw_i = k_silw * UPLIFT_i * PG_i * f_silw_i             # initial silicate weathering

# carbonate weathering, depend on runoff, vegetation, and co2, after COPSE reloaded
g_T_runoff_i = 1 + 0.087 * (temp_i - 288.15)             # runoff on carbonate
f_carb_i = g_T_runoff_i * ((1-VWmin_i) * k_planenhance * f_co2_abiotic_i + VEG_i * W_i)
carbw_i = k_carbw * PG_i * f_carb_i

# oxidw
oxw_fac_i = po2pal_i ** f_oxw_a
oxidw_i = k_oxidw * UPLIFT_i * oxw_fac_i 




# burial
mccb_i = silw_i + carbw_i                # dALK = 0
# locb_i = k_locb * VEG_i                  # continental burial
locb_i = oxidw_i * 0.5                   # continental and ocean burial in 1:1
mocb_i = oxidw_i * 0.5                   # dTC = 0

ap_i = oxidw_i  - locb_i - mocb_i + ccdeg_i - silw_i

# P cycle
p_i = p0 * np.sqrt(mocb_i/k_mocb)
# p_i = p0 * (mocb_i/k_mocb)
Pconc_i = (p_i/p0) * 2.2

# marine new production
newp_i = 117 * Pconc_i
# anoxic
ANOX_i =  1/(1+np.exp(-k_logistic * (k_uptake * (newp_i/newp0)-po2pal_i)))


# phosw_i = k_phosw * silw_i/k_silw 
mopb_i = mocb_i * ((ANOX_i/1000)+((1-ANOX_i)/250))

fepb_i = (k_fepb/k_oxfrac)*(1-ANOX_i)*(p_i/p0)

capb_i = k_capb * ((newp_i/newp0)**f_mocb_b)


# phosphorous balance
phosw_i = mopb_i+fepb_i+capb_i

# C isotope
phi_i = 0.01614 * (a_i / a0)   # fraction of C in atmosphere:ocean    


d_locb_i, D_P_i, d_mocb_i, D_B_i, d_mccb_i, d_ocean_i, d_atmos_i = fts.Cisotopefrac(temp_i, pco2pal_i, po2pal_i, phi_i)

delta_a_i = delta_mccb_i - d_mccb_i

d_mocb_i = -30
d_locb_i = -30

delta_g = -26
delta_c = 4
delta_vol = 0
# delta_mccb_i = 0
delta_g = (locb_i * ( delta_a_i+d_locb_i) + mocb_i * (delta_a_i+d_mocb_i) - ccdeg_i * delta_vol - carbw_i * delta_c + mccb_i * delta_mccb_i)/oxidw_i  
print(delta_g)
moldelta_a_i = delta_a_i * a_i

# U
u_i = u0 * 0.22
u_riv_i = u_riv0 * silw_i/k_silw
u_hydro_i = u_hydro0 * DEGASS_i 
u_anox_i = u_anox0 * (ANOX_i/0.0025) *u_i/u0
u_sed_i = u_sed0 * u_i/u0

u_riv_i - u_anox_i - u_sed_i - u_anox_i

# # assign initial values
# ystart = np.array([a_i, p_i, o_i, moldelta_a_i])
# t0 = -311e6
# tfinal = -285e6

# def derivs(t, y):
#     if t<-310e6:
#         t = -310e6
#     t_Ma = t/1e6
    
#     # read forcing
#     # forcing parameter
#     PG = ffunc.fPG(t_Ma)                   # Paleogeography effect on runoff/weathering, (Royer et al., 2014)
#     PG = ffunc.fPG(-310) 
#     DEGASS = ffunc.fhaq_D(t_Ma)            # metamorphic and volcanic degassing, inversion of sea-level curve. (Mills et al., 2017)
#     # DEGASS = ffunc.fhaq_D(-310)
#     UPLIFT = ffunc.fUPLIFT(t_Ma)           # tectonic uplift from sediment accumulation rates
#     # UPLIFT = ffunc.fUPLIFT(-310)                                              
#     EVO = ffunc.fEVO(t_Ma)                 # Land plant evolution
#     W = ffunc.fW(t_Ma)                     # Land plant enhancement of weathering
    
#     # retrieve the parameters
#     a, p, o, moldelta_a = y
    
#     delta_a = moldelta_a/a
#     # calcualte pco2
#     po2pal = o/o0
    
#     pco2pal = (a/a0)**2
#     pco2 = pco2atm0 * pco2pal
#     phi = 0.01614 * (a/a0)
    
#     temp =  k_CtoK + 15 + k_c * np.log(pco2pal) + k_l * t/570e6  + 1 
    
#     """--------------------- weathering  -------------------------""" 
#     # biomass
#     VEG = fts.land_biota(temp-k_CtoK, pco2, po2pal, EVO) 

#     # weathering rates, activation energy as 62 kJ/mol
#     f_T_silw = np.exp(0.09 * (temp-288.15))
#     f_T_ap = np.exp(0.0507 * (temp - 288.15))  # 35 kJ/mol

#     # runoff temperature dependence
#     f_T_runoff = ((1 + 0.038 * (temp - 288.15))**0.65)


#     f_co2_abiotic = pco2pal ** 0.5 

#     f_co2_biotic = (2*pco2pal / (1+pco2pal))**0.4

#     # o2 dependence
#     VWmin = min(VEG * W, 1)
#     f_silw = f_T_silw * f_T_runoff * ((1-VWmin) * k_planenhance * f_co2_abiotic + VEG * W)

#     silw = k_silw * UPLIFT * PG * f_silw             # initial silicate weathering

#     # carbonate
#     g_T_runoff = 1 + 0.087 * (temp - 288.15)             # runoff on carbonate
#     f_carb = g_T_runoff * ((1-VWmin) * k_planenhance * f_co2_abiotic + VEG * W)
#     carbw = k_carbw * PG * f_carb

#     # oxidw
#     oxw_fac = po2pal ** f_oxw_a
#     oxidw = k_oxidw * UPLIFT * oxw_fac

#     # degassing
#     ccdeg = k_ccdeg * DEGASS

#     # burial
#     mccb = silw + carbw
#     locb = k_locb * VEG
#     mocb = k_mocb * (p/p0) ** 2 

#     ap = ccdeg + oxidw  - locb - mocb - silw 
    
#     """--------------------- C isotope  -------------------------""" 
#     d_locb, D_P, d_mocb, D_B, d_mccb, d_ocean, d_atmos = fts.Cisotopefrac(temp, pco2pal, po2pal, phi)
    
#     delta_mccb = delta_a + d_mccb
#     moldelta_ap =  -locb*(delta_a+d_locb) - mocb * (delta_a+d_mocb) + oxidw*delta_g  + ccdeg*delta_vol + carbw*delta_c - mccb*delta_mccb 

#     """--------------------- P  -------------------------""" 
#     # P cycle
#     Pconc = (p/p0) * 2.2

#     # marine new production
#     newp = 117 * Pconc
#     # anoxic
#     ANOX =  1/(1+np.exp(-k_logistic * (k_uptake * (newp/newp0)-po2pal)))

#     # phosw_i = k_phosw * silw_i/k_silw 
#     mopb = mocb * ((ANOX/1000)+((1-ANOX)/250))

#     fepb = (k_fepb/k_oxfrac)*(1-ANOX)*(p/p0)
#     # fepb_i = (k_fepb/k_oxfrac)*(1-ANOX_i)
#     capb = k_capb * ((newp/newp0)**f_mocb_b)
#     # capb_i = k_capb * ((newp_i/newp0))

#     # phosphorous balance
#     phosw = phosw_i * (silw/silw_i)
    
#     pp = phosw - mopb-fepb-capb
    
#     yp = np.array([ap, pp, 0, moldelta_ap])
#     return yp

# """ --------------------------------------------------------  Main function ------------------------------------------------ """

# start_time = timeit.default_timer()

# print("\n@ Starting integration")
# print("[tstart tfinal]=[%.2e %.2e]\n" % (t0, tfinal))

# # ysol = solve_ivp(derives.derives,(t0,tfinal), ystart, t_eval= np.linspace(-313e6, -277e6, 400), method = 'LSODA')
# ysol = solve_ivp(derivs,(t0,tfinal), ystart,  method = 'LSODA')


# elapsed = timeit.default_timer() - start_time 

# print(f'{elapsed:.2f} used.')


# import pandas as pd
# df_sol = pd.DataFrame(ysol.y.T)
# df_sol.columns=['A',  'P',  'O', 'moldelta_A']


# df_sol['Age'] = ysol.t
# df_sol['phosphate_m'] = (df_sol['P']/k_oceanmass) * 1e6  # umol/kg
# df_sol['CO2_PAL'] = (df_sol['A']/a0)**2


# df_sol['CO2atm'] = df_sol['CO2_PAL'] * pco2atm0 * 1e6


# import matplotlib.pyplot as plt
# fig, axes = plt.subplots(figsize = (12,10), nrows = 3, ncols = 3)

# df_sol.plot(x='Age', y='CO2atm', ax=axes[0,0])
# plt.tight_layout()





