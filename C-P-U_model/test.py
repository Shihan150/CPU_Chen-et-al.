"""
A model for the C-P-U evolution 310-285 Ma
Author: Shihan Li
Several steps:
1. Auto spin to steady state in the beginning
2. Inverse to fit pco2 and d13c 
3. For d13c, randomly assume the other end-member d13c for co2 source
4. Output the U isotope value, compare with proxy records

Initial steady state: 
1. t = 310Ma
2. pCO2 = 500e-6
3. o = 4.4e19 (after COPSE results)
4. d13c = 4.42 (after proxy records)
5. d235u = -0.14 (after proxy records)

Forcing：
1. Paleogeography effect on runoff/weatheirng (Royer et al., 2014)
2. Metamorphic and volcanic degassing, inversion of sea-level curve. (Mills et al., 2017) ｜ or GEOCARBSURF (but low resolution)
3. tectonic uplift from sediment accumulation rates (Royer et al. 2014)
4. Land plant evolution on weathering. Always 1. Lenton et al., 2016
5. Plant enhancement of weathering. Always 1. Lenton et al., 2012

"""

# load packages
import numpy as np
from scipy.integrate import solve_ivp  
from scipy import interpolate
from scipy.optimize import fmin_l_bfgs_b, fmin
from tqdm import tqdm
import time

import timeit
import pandas as pd
import matplotlib.pyplot as plt


import load_forcings as ffunc         # COPSE_reloaded forcings
import functions as fts               # functions for ODE

'''--------------------------      load the proxy data  -------------------------------- '''
file_name = 'fit_data.xlsx'
df_u = pd.read_excel(file_name, sheet_name = 'U_AC_10%_LOESS')
df_d13c = pd.read_excel(file_name, sheet_name = 'C_AC_10%_LOESS')
df_pco2 = pd.read_excel(file_name, sheet_name = 'pCO2_LOESS')

df_pco2['age_new'] = -df_pco2['age']*1e6
df_d13c['age_new'] = -df_d13c['age']*1e6
df_u['age_new'] = -df_u['age']*1e6
df_u['d235u'] = df_u['mean']-0.27
# fig, axs = plt.subplots(2,2, figsize = (14,8))

# axs[0,0].plot(-df_u.age, df_u['mean'], marker = 'o', c='#8d2424', ls = '', label = 'd238U')
# axs[1,0].plot(-df_d13c.age, df_d13c['mean'], marker = 'o',c='#8d2424', ls='', label = 'd13c')
# axs[0,1].plot(-df_pco2['Age (Ma)'], df_pco2['LOESS CO2'], marker = 'o', c='#8d2424', ls='', label = 'pCO2')
# axs[1,1].remove()

# fig.supxlabel('Age (Ma)', fontsize=16, color="black")
# axs[0,0].set_ylabel('per mil', fontsize = 15)
# axs[1,0].set_ylabel('per mil', fontsize = 15)
# axs[0,1].set_ylabel('ppmv', fontsize = 15)

# axs[0,0].tick_params(axis="both", labelsize=10)
# axs[0,1].tick_params(axis="both", labelsize=10)
# axs[1,0].tick_params(axis="both", labelsize=10)

# axs[0,0].legend(fontsize = 15)
# axs[1,0].legend(fontsize = 15)
# axs[0,1].legend(fontsize = 15)

'''--------------------------      Model settings  -------------------------------- '''
# Some constants
# Allow some uncertainties for MC simulation later
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

# Modern Flux
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
   

# U, after Clarkson et al., 2018
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

'''--------------------------      Initial steady state  -------------------------------- '''
# initial values
t_i = -310e6
t_Ma_i = t_i/1e6

pco2_i = 500e-6                      # after proxy data
pco2pal_i = pco2_i/pco2atm0     
a_i = a0 * np.sqrt(pco2pal_i)
temp_i = k_CtoK + 15 + k_c * np.log(pco2pal_i) + k_l * t_i/570e6   # luminosity used imiplicitly
                                                                   # corrected 0 degrees by the sea surface temperature record
o_i = o0        # after COPSE result, equivalent to 25% O2atm, o2_pal 1.19

po2pal_i = o_i/o0    # relative o2 concentration 


# after proxy
delta_ocn_i = 4.42
delta_u_i = -0.14-0.27             # diagenetic correction, after Chen et al., 2018


#####     C cycle
# forcing parameter, degassing, uplift driving erosion, effects of vegetation and lithology on weatherability
PG_i = ffunc.fPG(t_Ma_i)                        # Paleogeography effect on runoff/weathering, (Royer et al., 2014)
DEGASS_i = ffunc.fhaq_D(t_Ma_i)                 # metamorphic and volcanic degassing, inversion of sea-level curve. (Mills et al., 2017) 
UPLIFT_i = ffunc.fUPLIFT(t_Ma_i)                # tectonic uplift from sediment accumulation rates
                                                # Royer et al. 2014
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

VWmin_i = min(VEG_i * W_i, 1)

f_silw_i = f_T_silw_i * f_T_runoff_i * ((1-VWmin_i) * k_planenhance * f_co2_abiotic_i + VEG_i * W_i)

# silw_i = k_silw * UPLIFT_i * PG_i * f_silw_i             # initial silicate weathering



# carbonate weathering, depend on runoff, vegetation, and co2, after COPSE reloaded
g_T_runoff_i = 1 + 0.087 * (temp_i - 288.15)             # runoff on carbonate
f_carb_i = g_T_runoff_i * ((1-VWmin_i) * k_planenhance * f_co2_abiotic_i + VEG_i * W_i)
carbw_i = k_carbw * PG_i * f_carb_i

# oxidw
oxw_fac_i = po2pal_i ** f_oxw_a
oxidw_i = k_oxidw * oxw_fac_i *DEGASS_i

# burial
mccb_i = silw_i + carbw_i                # dALK = 0
# locb_i = k_locb * VEG_i                  # continental burial
locb_i = oxidw_i * 0.5                   # continental and ocean burial in 1:1
mocb_i = oxidw_i * 0.5                   # dTC = 0

# ap_i = oxidw_i  - locb_i - mocb_i + ccdeg_i - silw_i

### P cycle, set weathering = burial
p_i = p0 * np.sqrt(mocb_i/k_mocb)
# p_i = p0 * (mocb_i/k_mocb)
Pconc_i = (p_i/p0) * 2.2

# marine new production
newp_i = 117 * Pconc_i
# anoxic
ANOX_i =  1/(1+np.exp(-k_logistic * (k_uptake * (newp_i/newp0)-po2pal_i)))



# phosw_i = k_phosw * silw_i/k_silw 
mopb_i = mocb_i * ((ANOX_i/1000)+((1-ANOX_i)/250))  # ocean burial

fepb_i = (k_fepb/k_oxfrac)*(1-ANOX_i)*(p_i/p0)      # fe-p 

capb_i = k_capb * ((newp_i/newp0)**f_mocb_b)        # ca3po4 burial


# phosphorous balance
phosw_i = mopb_i+fepb_i+capb_i


### C isotope, adjust the org weatheirng d13c
phi_i = 0.01614 * (a_i / a0)   # fraction of C in atmosphere:ocean    


d_locb_i, D_P_i, d_mocb_i, D_B_i, d_mccb_i, d_ocean_i, d_atmos_i = fts.Cisotopefrac(temp_i, pco2pal_i, po2pal_i, phi_i)

delta_a_i = delta_ocn_i - d_ocean_i
delta_mccb_i = delta_a_i + d_mccb_i

delta_mocb = -30
delta_locb = -30

delta_g = -26
delta_c = 4.5
delta_vol = 0
# delta_mccb_i = 0
# delta_g = (locb_i * ( delta_a_i+d_locb_i) + mocb_i * (delta_a_i+d_mocb_i) - ccdeg_i * delta_vol - carbw_i * delta_c + mccb_i * delta_mccb_i)/oxidw_i  
delta_g = (locb_i * ( delta_locb) + mocb_i * (delta_mocb) - ccdeg_i * delta_vol - carbw_i * delta_c + mccb_i * delta_mccb_i)/oxidw_i  

moldelta_a_i = delta_a_i * a_i



#### U, adjust the other
u_i = 1.96e13 
u_riv_i = u_riv0 * silw_i/k_silw
u_hydro_i = u_hydro0 * DEGASS_i 
u_anox_i = u_anox0 * (ANOX_i/0.0025) *u_i/u0
u_other_i = u_riv_i-u_hydro_i-u_anox_i

# isotope
delta_u_riv = -0.29
d_u_hydro = 0.2
delta_u_hydro = delta_u_i + d_u_hydro
d_u_anox = 0.6
delta_u_anox = delta_u_i + d_u_anox

delta_u_other = (u_riv_i * delta_u_riv - delta_u_hydro * u_hydro_i - delta_u_anox*u_anox_i)/u_other_i

d_u_other = delta_u_other-delta_u_i

moldelta_u_i = delta_u_i * u_i

ystart = np.array([a_i, p_i, o_i, moldelta_a_i, u_i, moldelta_u_i ])

def derivs(t, y, switch, pbar, state):
    
    t_original = t
    if t<-310e6:
        t = -310e6
        
    # t = -310e6
    t_Ma = t/1e6
    
    # read forcing
    # forcing parameter, allow for random sampling later
    PG = ffunc.fPG(t_Ma)/ffunc.fPG(-310)                   # Paleogeography effect on runoff/weathering, (Royer et al., 2014)
    PG = 1 
    DEGASS = ffunc.fhaq_D(t_Ma)/ffunc.fhaq_D(-310)              # metamorphic and volcanic degassing, inversion of sea-level curve. (Mills et al., 2017)
    DEGASS = 1
    UPLIFT = ffunc.fUPLIFT(t_Ma)/ffunc.fUPLIFT(-310)            # tectonic uplift from sediment accumulation rates
    UPLIFT = 1                                         
    EVO = ffunc.fEVO(t_Ma)/ffunc.fEVO(-310)                # Land plant evolution
    W = ffunc.fW(t_Ma)/ffunc.fEVO(-310)                     # Land plant enhancement of weathering
    
    
    
    # retrieve the parameters
    a, p, o, moldelta_a, u, moldelta_u = y
    
    delta_a = moldelta_a/a
    delta_u = moldelta_u/u
    # calcualte pco2
    po2pal = o/o0
    
    pco2pal = (a/a0)**2
    pco2 = pco2atm0 * pco2pal
    phi = 0.01614 * (a/a0)
    
    # temp = temp_i 
    temp = k_CtoK + 15 + k_c * np.log(pco2pal) + k_l * t/570e6  
    
    """--------------------- weathering  -------------------------""" 
    # biomass
    VEG = fts.land_biota(temp-k_CtoK, pco2, po2pal, EVO) 

    # weathering rates, activation energy as 62 kJ/mol
    f_T_silw = np.exp(0.09 * (temp-288.15))
    f_T_ap = np.exp(0.0507 * (temp - 288.15))  # 35 kJ/mol

    # runoff temperature dependence
    f_T_runoff = ((1 + 0.038 * (temp - 288.15))**0.65)


    f_co2_abiotic = pco2pal ** 0.5 

    f_co2_biotic = (2*pco2pal / (1+pco2pal))**0.4

    # o2 dependence
    VWmin = min(VEG * W, 1)
    f_silw = f_T_silw * f_T_runoff * ((1-VWmin) * k_planenhance * f_co2_abiotic + VEG * W)

    silw = silw_i * UPLIFT * PG * f_silw/f_silw_i            # initial silicate weathering

    # carbonate
    g_T_runoff = 1 + 0.087 * (temp - 288.15)             # runoff on carbonate
    f_carb = g_T_runoff * ((1-VWmin) * k_planenhance * f_co2_abiotic + VEG * W)
    carbw = carbw_i * PG * f_carb/f_carb_i

    # oxidw
    oxw_fac = po2pal ** f_oxw_a
    oxidw = oxidw_i * UPLIFT * oxw_fac/oxw_fac_i

    # degassing
    ccdeg = ccdeg_i * DEGASS

    # burial
    mccb = silw + carbw
    locb = locb_i * VEG/VEG_i
    mocb = mocb_i * (p/p_i) ** 2 

    ap = ccdeg + oxidw  - locb - mocb - silw + fcinp(t)*1e15/12
    
    """--------------------- C isotope  -------------------------""" 
    d_locb, D_P, d_mocb, D_B, d_mccb, d_ocean, d_atmos = fts.Cisotopefrac(temp, pco2pal, po2pal, phi)
    
    delta_mccb = delta_a + d_mccb
    delta_ocn = delta_a + d_ocean
    # moldelta_ap =  -locb*(delta_a+d_locb) - mocb * (delta_a+d_mocb) + oxidw*delta_g  + ccdeg*delta_vol + carbw*delta_c - mccb*delta_mccb 
    moldelta_ap =  -locb*(delta_locb) - mocb * (delta_mocb) + oxidw*delta_g  + ccdeg*delta_vol + carbw*delta_c - mccb*delta_mccb 


    """--------------------- P  -------------------------""" 
    # P cycle
    Pconc = (p/p0) * 2.2

    # marine new production
    newp = 117 * Pconc
    # anoxic
    ANOX =  1/(1+np.exp(-k_logistic * (k_uptake * (newp/newp0)-po2pal)))

    # phosw_i = k_phosw * silw_i/k_silw 
    mopb = mocb * ((ANOX/1000)+((1-ANOX)/250))

    fepb = (k_fepb/k_oxfrac)*(1-ANOX)*(p/p0)
    # fepb_i = (k_fepb/k_oxfrac)*(1-ANOX_i)
    capb = k_capb * ((newp/newp0)**f_mocb_b)
    # capb_i = k_capb * ((newp_i/newp0))

    # phosphorous balance
    phosw = phosw_i * (silw/silw_i)
    
    
    pp = phosw - mopb-fepb-capb
    
    # U cycle
    u_riv = u_riv_i * (silw/silw_i)
    u_hydro = u_hydro_i * DEGASS
    u_anox = u_anox_i * (ANOX/ANOX_i) * u/u_i
    u_other = u_other_i * (u/u_i)
    
    
    moldelta_up = u_riv * delta_u_riv - u_hydro*(delta_u+d_u_hydro) - u_anox*(delta_u+d_u_anox)- u_other*(delta_u+d_u_other)
    
    up = u_riv - u_hydro - u_anox - u_other
    
    
    
    
    if switch==1:
        yp = np.array([ap, pp, 0, moldelta_ap, up, moldelta_up])
        
        
        last_t, dt = state
        time.sleep(0.01)
        n = int((t_original-last_t)/dt)
        pbar.update(n)
        state[0] = last_t + dt*n
        
        
        return yp
    else:
        params = np.array([temp, ccdeg, oxidw, locb, mocb, silw, carbw, mccb, delta_ocn, phosw, mopb, fepb, capb, ANOX, UPLIFT, PG, DEGASS, VEG, EVO])
        return params

def cost_function(ems_new, tracer_type, *arg):
    
    if tracer_type == 'pco2' or tracer_type == 'pCO2':
        ems = np.array([arg[0], ems_new[0]])
        global fcinp
        fcinp = interpolate.interp1d([tems[0], tems[1]], [ems[0], ems[1]], bounds_error=False, fill_value=0)
        
        ysol = solve_ivp(derivs, (tems[0], tems[1]), y0, t_eval = np.array([tems[1]]), args={1},method = 'LSODA')
        
        a = ysol.y[0]
        pco2_model = (a/a0)**2 * pco2atm0 * 1e6
        errors = np.abs((pco2_model-f_target(tems[-1]))/f_target(tems[-1]))
        print(pco2_model)
        return errors

def run(Run_type, t0, tfinal, ystart, **kwargs):
    if Run_type == 1:
        # normal forward run
        global fcinp
        try:
            
            fcinp = interpolate.interp1d(kwargs['tems'],kwargs['ems'], fill_value = (0,0),bounds_error=False)  #
            
        except:
            fcinp = interpolate.interp1d((t0, tfinal),(0,0), fill_value = (0,0),bounds_error=False)  
        
        omega = 20
        
        with tqdm(total = 100, unit ='s') as pbar:
            
            ysol = solve_ivp(derivs,(t0,tfinal), ystart,  args=[1, pbar, [t0, (tfinal-t0)/100]], method = 'LSODA')
        
        # tidy the results
        t = ysol.t                # time
        y = ysol.y                # tracers
        
        nstep = len(t)

        params = np.zeros((nstep, 19))

        for i in range(nstep):
            params[i,:] = derivs(t[i], y[:,i], 0)
        
        df_params = pd.DataFrame(params)
        df_params.columns=['Temperature','ccdeg', 'oxidw', 'locb', 'mocb', 'silw', 'carbw', 'mccb', 'delta_ocn', 'phosw', 'mopb', 'fepb', 'capb', 'ANOX','UPLIFT', 'PG', 'DEGASS', 'VEG', 'EVO']
        df_params['Age'] = t
        
        df_sol = pd.DataFrame(ysol.y.T)
        df_sol.columns=['A',  'P',  'O', 'moldelta_A', 'U', 'moldelta_U']


        df_sol['Age'] = ysol.t
        df_sol['phosphate_m'] = (df_sol['P']/k_oceanmass) * 1e6  # umol/kg
        df_sol['U_m'] = (df_sol['U']/k_oceanmass)*1e6            # umol/kg
        df_sol['d235U'] = (df_sol['moldelta_U']/df_sol['U'])     # d235U
        df_sol['CO2_PAL'] = (df_sol['A']/a0)**2


        df_sol['CO2atm'] = df_sol['CO2_PAL'] * pco2atm0 * 1e6
        
        df_sol.to_csv("tracer.csv")
        df_params.to_csv("parameters.csv")
        
        return
    
    if Run_type == 2:
        # inverse calculation
        t_target = kwargs['target_t']
        tracer_target = kwargs['target_tracer']
        
        global f_target
        f_target = interpolate.interp1d(t_target,tracer_target, bounds_error=False) 
        
        tracer_eval = np.zeros(len(t_target))
        
        for i in range(len(t_target)-1):
            global y0
            if i:
                y0 = np.loadtxt('inv_ystart.dat')
            else:
                y0 = ystart
                
            global tems
        
            tems = t_target[i:i+2]
            ems = tracer_eval[i]
            
            try:
                init_guess = 2 * tracer_eval[i] - tracer_eval[i-1]
            except:
                init_guess = tracer_eval[i]
            
            ems_new = fmin(cost_function, init_guess, args = (kwargs['tracer_type'], ems), disp = True)
            
            tracer_eval[i+1] = ems_new
            
            # write out the ysol as the initial values for next loop
            ysol = solve_ivp(derivs, (tems[0], tems[1]), y0, t_eval = np.array([tems[1]]),args={1}, method = 'LSODA')
            np.savetxt('inv_ystart.dat', ysol.y)
        
        return np.vstack([t_target, tracer_eval])

# """ --------------------------------------------------------  Main function ------------------------------------------------ """

t0 = -311e6
tfinal = -285e6

# fitting target
t_target = -df_pco2['age'].to_numpy()*1e6
tracer_target = df_pco2['pco2'].to_numpy()


start_time = timeit.default_timer()


# cinp = run(2, t0, tfinal, ystart, tracer_type = 'pco2', target_t = t_target, target_tracer = tracer_target)

run(1, t0, tfinal, ystart, tems=cinp[0,:], ems=cinp[1,:])
print("\n@ Starting integration")
print("[tstart tfinal]=[%.2e %.2e]\n" % (t0, tfinal))

# run(1, t0, tfinal, ystart)


# print(f'{elapsed:.2f} used.')

# # 





# # ### plot
# df_sol = pd.read_csv("tracer.csv")
# df_params = pd.read_csv('parameters.csv')


# fig, axes = plt.subplots(figsize = (12,10), nrows = 2, ncols = 3)

# df_sol.plot(x='Age', y='CO2atm', ax=axes[0,0])
# df_pco2.plot(x='age_new', y='pco2', ax=axes[0,0], marker = '*', lw=0)
# df_sol.plot(x='Age', y='U_m', ax=axes[0,1])
# df_sol.plot(x='Age', y='d235U', ax=axes[0,2])
# df_u.plot(x='age_new', y = 'd235u', ax = axes[0,2], marker='*', lw=0)
# df_sol.plot(x='Age', y='phosphate_m', ax=axes[1,0])
# df_params.plot(x='Age', y ='ANOX', ax=axes[1,1])
# axes[1,2].remove()
# # df_params.plot(x='Age', y ='oxidw', ax=axes[1,2])
# # df_params.plot(x='Age', y ='locb', ax=axes[2,0])
# # df_params.plot(x='Age', y ='mocb', ax=axes[2,1])
# df_params.plot(x='Age', y ='silw', ax=axes[1,2])
# # df_params.plot(x='Age', y ='carbw', ax=axes[2,2])
# plt.tight_layout()

# fig, axes = plt.subplots(figsize = (12,10), nrows = 2, ncols = 3)

# df_params.plot(x='Age', y ='UPLIFT', ax=axes[0,0])
# df_params.plot(x='Age', y ='PG', ax=axes[0,1])
# df_params.plot(x='Age', y ='DEGASS', ax=axes[0,2])
# df_params.plot(x='Age', y ='VEG', ax=axes[1,0])
# df_params.plot(x='Age', y ='EVO', ax=axes[1,1])
# df_params.plot(x='Age', y ='ANOX', ax=axes[1,2])

# fig, ax = plt.subplots(figsize = (10,6))
# ax.plot(cinp[0,:],cinp[1,:])
# ax.set_ylabel('Gt C')
# ax.set_xlabel('Age')
