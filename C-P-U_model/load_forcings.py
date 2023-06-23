#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
load forcings
"""
import numpy as np
import pandas as pd
from scipy import interpolate

"""--------------------- forcing settings  -------------------------""" 

# 'PG' additional paleogeography forcing factor (extropolate)
PG_force = pd.read_excel('./forcings/royer_fDfAw.xlsx',  usecols = [0,1]).to_numpy()
fPG = interpolate.interp1d(-PG_force[:,0], PG_force[:,1],  fill_value = (1),bounds_error=False)

# DEGASS updated degassing based on sealevel inversion
haq_D_force = pd.read_excel('./forcings/D_haq_inversion_2017.xlsx',  usecols = [0,1]).to_numpy()
fhaq_D = interpolate.interp1d(-haq_D_force[:,0], haq_D_force[:,1],  fill_value = (1),bounds_error=False)

# UPLIFT updated uplift/erosion polynomial from geocarb 3 (extropolate)
UPLIFT_force = pd.read_excel('./forcings/berner_fr.xlsx', usecols = [0,1]).to_numpy()
fUPLIFT = interpolate.interp1d(-UPLIFT_force[:,0], UPLIFT_force[:,1],  fill_value = (1),bounds_error=False)


# CAL_NORM additional [Ca++] forcing factor, estimate = 'calnorm
CAL_NORM_force =  pd.read_excel('./forcings/Horita_Ca.xlsx',  usecols = [0,1]).to_numpy()
fCAL_NORM = interpolate.interp1d(np.hstack([-CAL_NORM_force[0,0]-1e-3, -CAL_NORM_force[:,0],-CAL_NORM_force[-1,0]+1e-3]), np.hstack([1, CAL_NORM_force[:,1],1]),fill_value = ([1]), bounds_error=False)

# new shale+coal+evaporite exposed area forcing which contributes to granite foring. estimate = 'orgevapnorm'
ORGEVAP_AREA_force =  pd.read_excel('./forcings/organics_and_evaporites_area.xlsx',  usecols = [0,1]).to_numpy()
fORGEVAP_AREA = interpolate.interp1d(np.hstack([-ORGEVAP_AREA_force[0,0]-1e-3, -ORGEVAP_AREA_force[:,0],-ORGEVAP_AREA_force[-1,0]+1e-3]), np.hstack([1, ORGEVAP_AREA_force[:,1],1]),fill_value = ([1]), bounds_error=False)

# new silicate exposed area forcing which contributes to granite forcing, estimate = silnorm
GRAN_force  = pd.read_excel('./forcings/shield_area.xlsx',  usecols = [0,2]).to_numpy()
fGRAN = interpolate.interp1d(np.hstack([-GRAN_force[0,0]-1e-3, -GRAN_force[:,0],-GRAN_force[-1,0]+1e-3]), np.hstack([1, GRAN_force[:,1],1]),fill_value = ([1]), bounds_error=False)

# COAL new coal depositional area forcing. estimate = 'coalnorm'
COAL_force = pd.read_excel('./forcings/coal_basin_frac_new.xlsx', usecols = [0,2]).to_numpy()
fCOAL = interpolate.interp1d(np.hstack([-COAL_force[0,0]-1e-3, -COAL_force[:,0],-COAL_force[-1,0]+1e-3]), np.hstack([1, COAL_force[:,1],1]),fill_value = ([1]), bounds_error=False)

# Newtiming of E
fEVO = interpolate.interp1d([-465,-445, -400, -350], [0,0.15,0.15,1],  fill_value = (0,1),bounds_error=False)

# New timing of W
fW = interpolate.interp1d([-465,-445, -400, -350], [0,0.75,0.75,1],  fill_value = (0,1),bounds_error=False)

# nEW FORCING OF c/p LAND TO COMBINE WITH BCOAL FORCING 
fCPland_relative = interpolate.interp1d([-465,-445, -345, -300], [1,2,2,1],  fill_value = (1),bounds_error=False)


# ramp up p weathering to get +2 per mil d13C plateau

fF_EPSILON = interpolate.interp1d([-465,-445, -410, -400], [1,1.5,1.5,1],  fill_value = (1),bounds_error=False)


# change B forcing to linear
fB = interpolate.interp1d([-150, 0], [0.75,1],  fill_value = (0.75,1),bounds_error=False)


# LIPS forcing
LIPs_version = 2  # LIPs_version: 1: mills2014g3; 2: copseRL; 3: copseRL_jw
if LIPs_version == 1:
    datafile = './forcings/ggge20620-sup-0001-suppinfol.xlsx'
    smoothempl = 0
    
elif LIPs_version == 2:
    datafile ='./forcings/CR_force_LIP_table.xlsx'
    smoothempl = 1
    
elif LIPs_version ==3:
    DATAFILE = './forcings/CR_force_LIP_table_jw.xlsx'
    smoothempl = 1

phanlip = pd.read_excel(datafile,  index_col = 0, skiprows=[1], nrows = 49,usecols=range(0,9))
phanlip.set_axis(['Names', 'Types', 'CFBareas', 'Allvolume', 'CO2min', 'CO2max', 'DegassDuration', 'PresentDayArea'], axis=1, inplace=True)



present_day_CFB_area = 4.8e6
present_day_OIB_area = 2e6
k_present_basalt_area = present_day_CFB_area + present_day_OIB_area

# forcing factor
RHO = 1              # additional enhancement to carbonate/silicate weathering
RHOSFW = 1           # additional enhancement to seafloor weathering
F_EPSILON = 1        # enhances nutrient weathering only
PG = 1               # PALEOGEOG WEATHERING FACTOR
BA = 1               # Basalt area relative to present day
GEOG = 1             # land temperature adjustment
GRAN_AREA = 1        # Granite area relative to present day
CARB_AREA = 1        # Granite area relative to present day
TOTAL_AREA = 1       # Total land area relative to present day
CPland_relative = 1  # C.P land multiplier
Ptoland = 1          # P to land multiplier
COAL = 1             # Coal depositional area multiplier
EVAP_AREA = 1        # Evaporite exposed area relative to present day
SALT = 1             # Evaporite depositional area multiplier
ORG_AREA = 1         # Shale + coal area relative to present day
SHALE_AREA = 1       # Shale area relative to present day
ORGEVAP_AREA = 1     # Shale + coal + evaporite area relative to present day


# CAL_NORM = 1       # prescribed calcium concentration - no default


co2pert = 0         # Pertubation (injection in mol.yr) to A reservoir 
co2pertmoldelta = 0

fireforcing = 1     # alteration of fire feedback

import numpy as np
import pandas as pd
from scipy.optimize import fsolve, root_scalar



# LIP forcing

class copse_load_phanlip:
    # Load LIPs excel file and constructs copse_force_LIP forcings for use in model
    
    def __init__(self, rev):
        self.rev = rev                    # data version: 1: mills2014g3, 2: copseRL, 3: copseRL_jw
        if rev == 1:
            self.datafile = './forcings/ggge20620-sup-0001-suppinfol.xlsx'
        elif rev == 2:
            self.datafile = './forcings/CR_force_LIP_table.xlsx'
        elif rev == 3:
            self.datafile = './forcings/CR_force_LIP_table_j2.xlsx'
        
        self.phanlip = pd.read_excel(self.datafile,  skiprows=[1], nrows = 49, usecols=range(0,9))
        self.phanlip.set_axis(['Allages', 'Names', 'Types', 'CFBareas', 'Allvolume', 'CO2min', 'CO2max', 'DegassDuration', 'PresentDayArea'], axis=1, inplace=True)
        
        
        
    
    def create_LIPs(self, co2releasefield, default_lambda, LIPtrange, smoothempl = 'false' ,areafac = 1, decayoffset = 0, default_co2releasetime = 2e5):
        
        self.phanlip['DegassDuration'] = self.phanlip['DegassDuration'].replace(np.nan, default_co2releasetime)
        self.phanlip = self.phanlip.fillna(0)
       
        
        LIPs = {}
        
        for i in range(len(self.phanlip['Allages'])):
            liptime = - self.phanlip['Allages'][i] * 1e6   # emplacement time (or NaN) in yr, relative to present-day = 0
            if co2releasefield == 'NoCO2':
                lipCO2 = 0
            else:
                lipCO2 = self.phanlib[co2releasefield][i]
            
            co2releasetime = self.phanlip['DegassDuration'][i] * 1e6 # myr -> yr
            
            peak_area = areafac * self.phanlip['CFBareas'][i]
            
            # calculate decay rate from present-day area, if available
            present_area = self.phanlip['PresentDayArea'][i]
            
            if present_area == 0:
                lamda = default_lambda
            
            else:
                lamda = np.log(present_area/peak_area)/(liptime + decayoffset)
            
            LIPs[i] = copse_force_LIP(liptime, self.phanlip['Names'][i], self.phanlip['Types'][i], \
                                      peak_area, smoothempl, lipCO2, decayoffset, lamda, co2releasetime)
        return LIPs
            
    def create_LIP_forcing(self, Tstart, Tend, co2releasefield, present_CFB_area, *argv):
        # create a lIP forcing function from data table, adjusting decay rates to match present_CFB_area
        #
        # Input
        # Tstart: start of time range where forcing required
        # Tend :  end of time range where forcing required
        # co2releasefield: 'NoCO2', 'CO2min', 'CO2max'
        # present_CFB_area: km2 target area, decay constant lambda (for those lips with unknown present day area) is adjusted to mett this
        # *argv: passed through to createdLIPs: LIPtrange, smoothempl, areafac, decayoffset, co2releasetime
        
        # Return
        # flips: composite forcing (basalt area and CO2), precalculated interpolation to grid for range Tstart:Tend
        # default_lambda: decay rate )for lips where present day area not available), to create present _CFB_area
        
        # solve for decay constant that provides present_CFB_area
        default_lambda = self.match_present_CFB_area(present_CFB_area, co2releasefield, *argv)
        # 1.1152839503606000e-8
        #self.match_present_CFB_area(present_CFB_area, co2releasefield, *argv)
        # print(default_lambda)
        # sys.exit()
        
        
        # create forcing
        # cell array of individual LIP forcings
        lipsAll = self.create_LIPs(co2releasefield, default_lambda, *argv)
        
        # create composite forcing, pre-calculated interpolation to tgrid
        # NB: forcing will only be interpolated and available for range Tstart - Tend
        tgrid = np.linspace(Tstart, Tend, num=int(1e4+1))
        flips = copse_force_composite(lipsAll, tgrid, ['CFB_area', 'LIP_CO2', 'LIP_CO2moldelta'], [0, 0 ,0])
        flips.fastinterp = 1
        
        # # check: evaluate present-day CFB area
        # D = {}
        # D['CFB_area'] = 0
        # D['LIP_CO2'] = 0
        # D['LIP_CO2moldelta'] = 0
        
        # D = flips.force(init.time_present_yr, D)
        
        # frac_area_err = abs(D['CFB_area'] - present_CFB_area)/present_CFB_area
        # if frac_area_err > 1e-8:
        #     print('LIP present_day_area failed')
        #     sys.exit()
        
        return  default_lambda, lipsAll, flips
            
    def match_present_CFB_area(self, present_CFB_area, co2releasefield, *argv):
        # solve for decay rate that produces specified present_CFB_area (km2)
        # *argv is passed through to createLIPs: smoothempl, areafac, decayoffset, co2releasetime
        areaDiff = lambda lamb: self.get_present_CFB_area(co2releasefield, lamb, *argv) - present_CFB_area
        default_lambda = fsolve(areaDiff, (0.25e-8))
        # default_lambda = root_scalar(areaDiff, method = 'bisect',  bracket = [0.25e-8, 4e-8])
        return default_lambda

    
    def get_present_CFB_area(self, co2releasefield, default_lambda, *argv):
        lipsAll = self.create_LIPs(co2releasefield, default_lambda, *argv)
        
        
        flips = copse_force_composite(lipsAll)
        
        D = {}
        D['CFB_area'] = 0
        D['LIP_CO2'] = 0
        D['LIP_CO2moldelta'] = 0
        
        D = flips.force(0, D)
        
      
        
        return D['CFB_area']
    
    

class copse_force_composite :
    tgrid = []
    
    def __init__(self, forcings,  *argv, fastinterp = 0):
        self.forcings = forcings 
        if len(argv):
            self.tgrid = argv[0]
            self.dfields = argv[1]
            self.dfieldstoset = argv[2]
            self.calcgrid()
            
        
    def force(self, tforce_or_tmodel, D):
        if self.tgrid:
            D = self.force_grid(tforce_or_tmodel, D)
        else:
            
            D = self.force_individual(tforce_or_tmodel, D)
        
        return D
    
    def force_individual(self, tforce_or_tmodel, D):
        
        for i in range(len(self.forcings)):
            
            D = self.forcings[i].force(tforce_or_tmodel, D)
            
        return D
    
    def calcgrid(self):
        D = {}
        
        for i in range(len(self.dfields)):
            D[self.dfields[i]] = np.zeros(len(self.tgrid))            
        self.Dgrid = self.force_individual(self.tgrid, D)
    
    # def force_grid(self, tforce_or_tmodel, D):
    #     for i in length(self.dfields):
    #         if self.fastinterp:
    #             interpval
    #     return D


        
 
class copse_force_LIP:
    
    smoothempl = ''
    
    def __init__(self, liptime, Name, Type, peakCFBarea, smoothempl, co2_potential, decayoffset, lamb, co2releasetime):
        self.liptime = liptime                   # yr
        self.Name = Name                         # vector of spreadshheet rows
        self.Type = Type                         
        self.peakCFBarea = peakCFBarea           # km2 peak area exposed to subaerial weathering       
        self.smoothempl = smoothempl             # km2 present-day area
        self.co2_potential = co2_potential       # co2 released (mol C)
        self.co2_d13c = -5                       # co2 release isotopic composition todo value
        
        self.decayoffset = decayoffset           # yr  wait before erosion begins
        self.lamb = lamb                         # 1/yr lambda for exponential area decay
        self.co2releasetime = co2releasetime    # yr timeframe for co2 release (yrs) (full width of gaussian) 
        
    
    def calc_LIP_CFBarea(self, timeyr):
        smoothing = 3e-6
        offset = 1.5e6
        sigmoid = lambda x,y : 1/(1 + np.exp(-y * x))
        
        if self.smoothempl:
            
            empl = sigmoid(smoothing, timeyr - self.liptime + offset)
        
        else:
            empl = 0.5 + 0.5 * np.sign(timeyr - self.liptime)
        
        # exponential decay
        erode = (0.5 + 0.5 * np.sign(timeyr - self.liptime - self.decayoffset)) * (1-np.exp(-self.lamb * (timeyr - self.liptime - self.decayoffset)))
        
        area = self.peakCFBarea * (empl - erode)
        
        return area
    
    def calc_LIP_co2(self, timeyr):
        gaussnorm = lambda x, y: 1/(x * (2 * np.pi) ** 0.5) * np.exp(-y ** 2/ ( 2 * x ** 2))
        co2release = self.co2_potential * gaussnorm(0.5 * self.co2releasetime, timeyr - self.liptime)
        
        return co2release
    
    def force(self, tforce_presentdayiszeroyr, D):
        # calculate basalt area and co2 relate
        D['CFB_area'] += self.calc_LIP_CFBarea(tforce_presentdayiszeroyr)
        
        co2 = self.calc_LIP_co2(tforce_presentdayiszeroyr)
        D['LIP_CO2'] += co2
        # D['LIP_CO2'] + co2
        D['LIP_CO2moldelta'] += co2 * self.co2_d13c
        
        return D
        
        
LIPs_version = 2  # LIPs_version: 1: mills2014g3; 2: copseRL; 3: copseRL_jw
if LIPs_version == 1:
    datafile = './forcings/ggge20620-sup-0001-suppinfol.xlsx'
    smoothempl = 0
    
elif LIPs_version == 2:
    datafile ='./forcings/CR_force_LIP_table.xlsx'
    smoothempl = 0
    
elif LIPs_version ==3:
    DATAFILE = './forcings/CR_force_LIP_table_jw.xlsx'
    smoothempl = 1


present_day_CFB_area = 4.8e6    
phanlip = copse_load_phanlip(2)
default_lambda, lipsAll, flips = phanlip.create_LIP_forcing(-1e9, 0, 'NoCO2', present_day_CFB_area,[-np.inf, np.inf], smoothempl)

def fLIP(lipsAll , t):
    flips = copse_force_composite(lipsAll)
    D = {}
    D['CFB_area'] = 0
    D['LIP_CO2'] = 0
    D['LIP_CO2moldelta'] = 0
    
    D = flips.force(t, D)
    
  
    
    return D['CFB_area'], D['LIP_CO2'], D['LIP_CO2moldelta']


