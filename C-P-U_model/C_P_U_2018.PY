"""
A model for the C-P-U evolution 310-290 Ma
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
1. linear weatherability scale for silicate weathering
2. linear weatherability scale for carbonate weathering
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy import interpolate
from scipy.optimize import fmin_l_bfgs_b, fmin

import timeit
import pandas as pd
import matplotlib.pyplot as plt

import functions as fts               # functions for ODE
import emcee


### load the target data
target = pd.read_csv('target.csv')

'''--------------------------    Initial fluxes at t = 310ma -------------------------------- '''

# Constants
k_logistic = 12              # new determines slope of logistic anoxia function, COPSE_reloaded, Clarkson et al. 2018
k_uptake = 0.5000            # new determines efficiency of nutrient uptake in anoxia function; COPSE_reloaded; Clarkson et al. 2018
k_CtoK = 273.15    # convert degrees C to Kelvin
k_c = 4.3280                 # new determines climate sensitivity to CO2
k_l = 7.4000                 # new determines temperature sensitivity to luminosity
k_oxfrac = 0.9975            # updated initial oxic fraction
k_oceanmass = 1.397e21       # ocean mass (kg)

f_oxw_a = 0.5                # oxidative weathering dependency on O2 concentration
f_mocb_b = 2                 # marine organic carbon burial dependency on new production

a0 = 3.193e18       # atmosphere-ocean co2
o0 = 3.7e19         # atmosphere-ocean o2

pco2_i = 500e-6      
pco2atm0 = 280e-6
pco2pal_i = pco2_i/280e-6
a_i = a0 * np.sqrt(pco2pal_i)

delta_ocn_i = 4.46
delta_u_i = -0.16-0.27             # diagenetic correction, after Chen et al., 2022

k_oxidw = 5e12      # oxidative weatheirng + degassing, kept constant for simplicity
k_locb = 2.5e12     # continental organic C burial
k_mocb = 2.5e12     # marine organic C burial

# carbon isotope
delta_mocb = -30
delta_locb = -30


delta_c = 4.5
delta_vol = -1
# organic weathering is kept free to close the 13c cycle

# organic carbon cycle follows the modern world
oxidw_i = k_oxidw
locb_i = k_locb
mocb_i = k_mocb

k6_fepb = 1e10               # updated Fe-P burial (mol/year)
k7_capb = 2e10               # updated Ca-P burial (mol/year)
k_mopb = 1e10                # organic-P burial in the ocean
k10_phosw = 4e10             # updated P weathering (mol/year)

newp0 = 117 * 2.2            # new production (umol/kg)
p0 = 2.2 * 1e-6 * k_oceanmass               # ocean (phosphate) phosphorus




# U cycle follows the modern value, after clarkson et al., 2018
u0 = 1.85e13         # modern U in the ocean

u_riv0 = 4.79e7       # river input
u_hydro0 = 5.7e6      # hydrothermal output
u_anox0 = 6.2e6       # anoxic sink
u_other0 = 3.6e7

u_i = u0
u_riv_i = u_riv0
u_hydro_i = u_hydro0
u_anox_i = u_anox0
u_other_i = u_other0

delta_u_riv = -0.29
d_u_hydro = 0.2
delta_u_hydro = delta_u_i + d_u_hydro
d_u_anox = 0.6
delta_u_anox = delta_u_i + d_u_anox

delta_u_other = (u_riv_i * delta_u_riv - delta_u_hydro * u_hydro_i - delta_u_anox*u_anox_i)/u_other_i

d_u_other = delta_u_other-delta_u_i

# isotope fractionations
# d238u_riv = -0.29
# u_frac_anox = 0.5
# u_frac_sed = 0.0156
# u_frac_hydro = 0.2



###########       Parameters for inversion       ################
silw_i = 5              # e12 mol/yr, 5-10.5
carbw_i = 7            # e12 mol/yr, 4-14

temp_i = 2.85         # K, 283-287

po2pal_i = 1         # 1.0 - 1.2

ppal_i = 1.4         # 1.0 - 1.5

scale_silw = [1.4, 1.6]  # integrated long-term silicate weathering scale at t =  305, 300, 295, 290 Ma
scale_carbw = 1.2 # integrated long-term carbonate weathering modifier at t = 305, 300, 295, 290 Ma
scale_degassing = 1.0   # relative outgassing scale, 0.8-1.2
scale_u_riv = 0.2   # relative d_u_riv modifier to fit the lont-term d238u trend
scale_d13c_oxidw = 0.7  # relative oxidw modifier to fit the long-term d13c trend

alpha = 0.33             # co2 dependence, 0.2-0.5
te = 3.4                 # e-folding temperature dependence of continental weathering, 5-50, after Krissansen-Totton et a., 2018

p_on_silw = 1          # P dependence on silw flux

cinput = np.array([1.5,2,1.8,4])   # *1e15 Carbon for 4 carbon emission events, corresponding time is manually defined

params =[silw_i, carbw_i, temp_i, po2pal_i, ppal_i, scale_silw[0], scale_silw[1], scale_carbw,
         # scale_carbw[0],scale_carbw[1]
         scale_degassing, scale_u_riv, scale_d13c_oxidw, alpha, te, 
         # p_on_silw,
         cinput[0], cinput[1], cinput[2], cinput[3]]
# np.array([silw_i, carbw_i, temp_i, po2pal_i, scale_silw, scale_carbw, scale_degassing, te, cinput])



'''--------------------------    Probabilty function   -------------------------------- '''
def log_probability(params):
    lp = log_prior(params)
    
    if not np.isfinite(lp):
        return -np.inf
    
    else:
        print(log_prior(params), params)
        global fscale_silw, fscale_carbw, temp_i, fscale_degassing,silw_i, carbw_i, alpha, te, ccdeg_i,o_i,p_i, fcinp, ANOX_i, phi_i, ppal_i, fscale_u_riv, fscale_d13c_oxidw, p_on_silw
        
        cinput = np.zeros(4)
        scale_silw = np.zeros(2)
        silw_i, carbw_i, temp_i, po2pal_i, ppal_i, scale_silw[0], scale_silw[1], scale_carbw, scale_degassing, scale_u_riv, scale_d13c_oxidw, alpha, te,  cinput[0], cinput[1], cinput[2], cinput[3]= params
        p_on_silw = 1.2
        silw_i*=1e12
        carbw_i*=1e12
        temp_i = temp_i*100
        te *= 10
        cinput *= 10
        
    
        
       
        # interp_scale_carbw.extend(scale_carbw)
        
        
        
        fscale_silw = interpolate.interp1d([-310e6, -300e6, -290e6],[1,scale_silw[0], scale_silw[1]], bounds_error = False, fill_value = 1)
        fscale_carbw = interpolate.interp1d([-310e6,-290e6], [1,scale_carbw], bounds_error = False, fill_value = 1)
        fscale_degassing = interpolate.interp1d([-310e6, -290e6],[1,scale_degassing], bounds_error = False, fill_value = 1)
        fscale_u_riv = interpolate.interp1d([-310e6, -290e6],[1,scale_u_riv], bounds_error = False, fill_value = 1)
        fscale_d13c_oxidw = interpolate.interp1d([-310e6, -290e6],[1,scale_d13c_oxidw], bounds_error = False, fill_value = 1)
    
        cinput_age = [-305.3e6, -303.8e6, -302.30e6, -301.30e6, -298.32e6, -296.66e6, -295.13e6, -293.73e6]
        # cinput_age = [-304.3e6, -304.15e6, -299.19e6, -298.43e6, -297e6, -296.4e6, -295.73e6, -293.91e6]
        cinput_rate = np.array(cinput)* 1e18 /np.array([-303.84e6 + 304.3e6, -301.30e6+302.30e6, -296.4e6+297e6, 295.73e6-293.91e6])    
        # cinput_rate = np.array(cinput)* 1e18 /np.array([-304.15e6 + 304.3e6, -298.43e6+299.19e6, -296.4e6+297e6, 295.73e6-293.91e6])    
        
        fcinp = interpolate.interp1d(cinput_age, [cinput_rate[0], 0, cinput_rate[1], 0, cinput_rate[2],0, cinput_rate[3],0], kind = 'zero', bounds_error = False, fill_value = 0)
        
        t_eval = np.sort(target.age.values)
        np.savetxt('test.dat', fcinp(t_eval)/12/1e13)
        
        ##############   initialize the model at t = 310Ma #################
        
        # Carbon
        mccb_i = silw_i + carbw_i    # total carbon burial to maintain the alkalinity balance, after COPSE
        o_i = o0 * po2pal_i
    
        
        ccdeg_i = silw_i
        
        # ap_i = ccdeg_i+oxidw_i-locb_i-mocb_i+carbw_i-mccb_i       # check the balance of carbon cycle
        
        # C isotope
        phi_i = 0.01614 * (a_i/a0)  # fraction of C in atmosphere:ocean
        
        d_locb_i, D_P_i, d_mocb_i, D_B_i, d_mccb_i, d_ocean_i, d_atmos_i = fts.Cisotopefrac(temp_i, pco2pal_i, po2pal_i, phi_i)
        
        delta_a_i = delta_ocn_i - d_ocean_i
        delta_mccb_i = delta_a_i + d_mccb_i
        
        global delta_g
        delta_g = (locb_i * ( delta_locb) + mocb_i * (delta_mocb) - ccdeg_i * delta_vol - carbw_i * delta_c + mccb_i * delta_mccb_i)/oxidw_i
        
        moldelta_a_i = a_i * delta_a_i
        
        # moldelta_ap_i = -ccdeg_i*delta_vol - oxidw_i*delta_g + locb_i*delta_locb + mocb_i*delta_mocb - carbw_i*delta_c + mccb_i*delta_mccb_i 
        
        
        # P cycle
        global p_i, newp_i, ANOW_i, mopb_i, fepb_i, capb_i, phosw_i
        p_i = p0 * ppal_i
        newp_i = 117 * (p_i/p0) * 2.2
        ANOX_i = 1/(1+np.exp(-k_logistic * (k_uptake * (newp_i/newp0)-po2pal_i)))   
        mopb_i = mocb_i * ((ANOX_i/1000)+((1-ANOX_i)/250))  # ocean burial
        fepb_i = (k6_fepb/k_oxfrac)*(1-ANOX_i)*(p_i/p0)
        capb_i = k7_capb * ((newp_i/newp0)**f_mocb_b)
        
        phosw_i = mopb_i + fepb_i + capb_i
        
        pp_i = phosw_i - mopb_i -fepb_i - capb_i
        
        # U cycle
        # up = u_riv_i - u_hydro_i - u_anox_i - u_other_i
        # moldelta_up_i = u_riv_i*d238u_riv - u_hydro_i*delta_u_hydro - u_anox_i*delta_u_anox - u_other_i*delta_u_other
        # print(up)
        # print(moldelta_up_i)
        moldelta_u_i = u_i * delta_u_i
        
        ystart = np.array([a_i, p_i, o_i, moldelta_a_i, u_i, moldelta_u_i])
        t0 = -310e6
        tfinal = -290e6
    
    
    
        start_time = timeit.default_timer()
        
        ysol = solve_ivp(derivs,(t0,tfinal), ystart, args={1}, method = 'LSODA', t_eval = t_eval, max_step = 1e4)
        # ysol = derivs(-310e6, ystart, 1)
        
        # print("\n@ Starting integration")
        # print("[tstart tfinal]=[%.2e %.2e]\n" % (t0, tfinal))
        
        if np.isnan(ysol.y).any():
            return -np.inf
        
        else:
            t = ysol.t                # time
            y = ysol.y                # tracers

            nstep = len(t)

            params = np.zeros((nstep, 14))

            for i in range(nstep):
                params[i,:] = derivs(t[i], y[:,i], 0)

            df_params = pd.DataFrame(params)
            df_params.columns=['Temperature','ccdeg', 'oxidw', 'locb', 'mocb', 'silw', 'carbw', 'mccb', 'delta_ocn', 'phosw', 'mopb', 'fepb', 'capb', 'ANOX']
            df_params['Age'] = t

            df_sol = pd.DataFrame(ysol.y.T)
            df_sol.columns=['A',  'P',  'O', 'moldelta_A', 'U', 'moldelta_U']


            df_sol['Age'] = ysol.t
            df_sol['phosphate_m'] = (df_sol['P']/k_oceanmass) * 1e6  # umol/kg
            df_sol['p/p0'] = df_sol['phosphate_m']/2.2
            df_sol['U_m'] = (df_sol['U']/k_oceanmass)*1e6            # umol/kg
            df_sol['d235U'] = (df_sol['moldelta_U']/df_sol['U'])     # d235U
            df_sol['CO2_PAL'] = (df_sol['A']/a_i)**2
            df_sol['d13c'] = (df_sol['moldelta_A']/df_sol['A'])      # d13c

            df_sol['CO2atm'] = df_sol['CO2_PAL'] * pco2_i * 1e6
            df_sol['O2PAL'] = df_sol['O']/o0
            # df_sol.to_csv("tracer.csv")
            # df_params.to_csv("parameters.csv")

            fig, axes = plt.subplots(figsize = (12,10), nrows = 3, ncols = 3)

            df_sol.plot(x='Age', y='CO2atm', ax=axes[0,0])
            target.plot(x='age', y='pco2', ax=axes[0,0], marker = '*', lw=0)
            df_sol.plot(x='Age', y='U_m', ax=axes[0,1])
            df_sol.plot(x='Age', y='d235U', ax=axes[0,2])
            target.plot(x='age', y = 'u', ax = axes[0,2], marker='*', lw=1.5)
            df_sol.plot(x='Age', y='phosphate_m', ax=axes[1,0])
            df_params.plot(x='Age', y ='ANOX', ax=axes[1,1])
            # axes[1,2].remove()
            # df_params.plot(x='Age', y ='oxidw', ax=axes[1,2])
            # df_params.plot(x='Age', y ='locb', ax=axes[2,0])
            # df_params.plot(x='Age', y ='mocb', ax=axes[2,1])
            df_params.plot(x='Age', y ='silw', ax=axes[1,2])
            # df_params.plot(x='Age', y ='carbw', ax=axes[2,2])
            df_sol.plot(x='Age', y='d13c', ax=axes[2,0])
            target.plot(x='age', y='d13c', ax=axes[2,0], marker = '*', lw=1.5)
            df_sol.plot(x='Age', y='O2PAL', ax=axes[2,1])
            
            axes[2,2].plot(t, fcinp(t)/1e15)
            axes[2,2].set_ylabel('Carbon emission (Gt/yr)')


            plt.tight_layout()

            pco2_proxy = target['pco2'].values
            pco2_std = target['pco2_std'].values
            d13c_proxy = target['d13c'].values
            d13c_std = target['d13c_std'].values
            u_proxy = target['u'].values
            u_std = target['u_std'].values

            pco2_model = df_sol['CO2atm'].values
            d13c_model = df_sol['d13c'].values
            u_model = df_sol['d235U'].values

            u_index1 = range(21)
            u_index2 = 35
            
            index = [0,5,11,15,20,37, 44,50,55,60,65,70,73]

            sum_diff = sum((pco2_proxy[::3]-pco2_model[::3])**2/(pco2_std[::3]**2)) +   sum((u_proxy[index]-u_model[index])**2/(u_std[index]**2)) 
    #         sum((u_proxy-u_model)**2/(u_std**2))
            # 0.2 *  sum((d13c_proxy-d13c_model)**2/(d13c_std**2)) +  0.5
            

           
            print(-0.5*sum_diff)

            return lp-0.5 * sum_diff



def derivs(t, y, switch):
    # if t<-310e6:
    #     t = -310e6

    # t = -310e6
    t_Ma = t/1e6

  


    # retrieve the parameters
    a, p, o, moldelta_a, u, moldelta_u = y

    delta_a = moldelta_a/a
    delta_u = moldelta_u/u
    # calcualte pco2
    po2pal = o/o0

    pco2pal = (a/a_i)**2
    pco2 = pco2_i * pco2pal
    phi = 0.01614 * (a/a_i)

    # temp = temp_i
    temp = temp_i + k_c * np.log(pco2pal) + k_l/570e6 * (t+310e6)
   
    
    diff_temp = temp - temp_i

    """---------------------  Carbon  -------------------------"""
    
    # degassing
    ccdeg = ccdeg_i * fscale_degassing(t)
   
    silw = fscale_silw(t) * silw_i * (pco2pal ** alpha) * np.exp(diff_temp/te) 
    
    
    
    
    carbw = fscale_carbw(t) * fscale_silw(t) * carbw_i * (pco2pal ** alpha) * np.exp(diff_temp/te)          

    # oxidw
    oxw_fac = po2pal ** f_oxw_a
    oxidw = oxidw_i *  oxw_fac 

    

    # burial
    mccb = silw + carbw
    locb = locb_i * (2*pco2pal/(1+pco2pal))
    mocb = mocb_i * (p/p_i) ** 2

    ap = ccdeg + oxidw  - locb - mocb - silw + fcinp(t)/12
   
   

    """--------------------- C isotope  -------------------------"""
    d_locb, D_P, d_mocb, D_B, d_mccb, d_ocean, d_atmos = fts.Cisotopefrac(temp, pco2pal * pco2pal_i, po2pal*po2pal_i, phi * np.sqrt(pco2pal_i))

    delta_mccb = delta_a + d_mccb
    delta_ocn = delta_a + d_ocean
    # moldelta_ap =  -locb*(delta_a+d_locb) - mocb * (delta_a+d_mocb) + oxidw*delta_g  + ccdeg*delta_vol + carbw*delta_c - mccb*delta_mccb
    moldelta_ap =  -locb*(delta_locb) - mocb * (delta_mocb) + oxidw*delta_g*fscale_d13c_oxidw(t) + ccdeg*delta_vol + carbw*delta_c - mccb*delta_mccb + fcinp(t)/12 * -20
    
    
    

    """--------------------- P  -------------------------"""
    # P cycle
    Pconc = (p/p0) * 2.2

    # marine new production
    newp = 117 * Pconc
    # anoxic
    ANOX =  1/(1+np.exp(-k_logistic * (k_uptake * (newp/newp0)-po2pal)))
    
    # phosw_i = k_phosw * silw_i/k_silw
    mopb = mocb * ((ANOX/1000)+((1-ANOX)/250))

    fepb = (k6_fepb/k_oxfrac)*(1-ANOX)*(p/p0)
    # fepb_i = (k_fepb/k_oxfrac)*(1-ANOX_i)
    capb = k7_capb * ((newp/newp0)**f_mocb_b)
    # capb_i = k_capb * ((newp_i/newp0))

    # phosphorous balance
    phosw = phosw_i * ((silw)/(silw_i)) ** p_on_silw


    pp = phosw-mopb-fepb-capb
    
    """--------------------- O  -------------------------"""
    op = locb + mocb - oxidw 
    
    # U cycle
    u_riv = u_riv_i * (silw/silw_i) 
    u_hydro = fscale_degassing(t) *u_hydro_i 
    u_anox = u_anox_i * (ANOX/ANOX_i) * u/u_i
    
    
    
    u_other =  u_other_i * (u/u_i)
    


    moldelta_up = u_riv * delta_u_riv * fscale_u_riv(t) - u_hydro*(delta_u+d_u_hydro) - u_anox*(delta_u+d_u_anox)- u_other*(delta_u+d_u_other)

    up = u_riv - u_hydro - u_anox - u_other



    
    if switch:
        yp = np.array([ap, pp, op, moldelta_ap, up, moldelta_up])
        return yp
    else:
        params = np.array([temp, ccdeg, oxidw, locb, mocb, silw, carbw, mccb, delta_ocn, phosw, mopb, fepb, capb, ANOX])
        return params

def log_prior(theta):
    
    cinput = np.zeros(4)
    scale_silw = np.zeros(2)
    silw_i, carbw_i, temp_i, po2pal_i, ppal_i, scale_silw[0], scale_silw[1], scale_carbw, scale_degassing, scale_u_riv, scale_d13c_oxidw, alpha, te,cinput[0], cinput[1], cinput[2], cinput[3]= theta
    if 2<=silw_i<=12 and 4<=carbw_i<=14 and 2.83<=temp_i<=2.87 and 1.0<=po2pal_i<=1.3 and 1.0<=ppal_i<=1.5 and 1.0<=scale_silw.any()<=2.0 and 1.0<= scale_carbw<=2 and 0.8<=scale_degassing<=1.2 and 0<=scale_u_riv<=1.0 and 0.2 <= scale_d13c_oxidw<=1.0 and 0.2<=alpha<=0.5 and 0.5<=te<=5 and 0.5<=cinput[0]<=3 and 0.5<=cinput[1]<=3 and 0.5<=cinput[2]<=3 and .5<=cinput[3]<=4  :
        return 0.0
    return -np.inf

    
   
log_probability([5.09610023, 9.079668,   2.84611144, 1.00953257, 1.42505597, 1.4, 2.0,
1, 1.07570367, 0.64767465, 0.82393394, 0.36004522, 1.88942291, 1.5,
1.5, 1, 3.47145794])
