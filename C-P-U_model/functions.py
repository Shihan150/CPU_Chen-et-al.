#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for the C-P-U model
"""
import numpy as np

def land_biota(temp, pco2, po2pal, EVO):
    # COPSE_LANDBIOTA COPSE OCT dynamic land vegetation model
    # represents land plant coverage and biomass, and organic carbon thus produced.
    # effect of temp on VEG, land primary productivity
    # depend on pCO2, pO2, and global temperature
    # fire feedback not considered here
    
    # Input: global temp, pco2, po2, EVO scalar
    # Output: VEG
    
    V_T = 1 - ((temp - 25)/25) ** 2  # maximum productivity at 25 C, decreasing to 0 at 0 and 50 C
    
    
    # effect of CO2 on VEG, Michaelis-Menten equation (Volk, 1987)
    P_atm = pco2 * 1e6 
    P_half = 183.6          # p_half, adopted so that the maximum productivity possible is twice pre-industrial levels                          
    P_min = 10              # p_min is the minimum pCO2 under which plants can grow
    V_co2 = (P_atm - P_min) / (P_half + P_atm - P_min)
    
    
    
    # EFFECT OF O2 on VEG, Refield Revisited (LW2) model
    V_o2 = 1.5 - 0.5 * po2pal
    
    # full VEG limitation, original 
    V_npp = 2 * EVO * V_T * V_o2 * V_co2
    
    
    return V_npp

def Cisotopefrac(Tkelvin, pCO2PAL, pO2PAL, phi = 0.01614):
    # input
    # phi: fraction of C in atmosphere:ocean
    
    # ocean total dissolved co2
    d_ocean = phi * (9483/Tkelvin - 23.89)
    # atmosphere CO2
    d_atmos = (phi-1) * (9483/Tkelvin-23.89)
    
    # marine calcite burial
    d_mccb = d_ocean + 15.10 - 4232/Tkelvin
    
    # fractionation between marine organic and calcite burial
    D_B = 33 - 9/np.sqrt(pCO2PAL) + 5 * (pO2PAL-1)
    # marine organic carbon burial
    d_mocb = d_mccb - D_B
    
    # fractionation between terrestrial organic burial and atmospheric CO2
    D_P = 19 + 5 * (pO2PAL - 1)
    d_locb = d_atmos - D_P
    
    return d_locb, D_P, d_mocb, D_B, d_mccb, d_ocean, d_atmos
