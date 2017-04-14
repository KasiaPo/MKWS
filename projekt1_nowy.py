# -*- coding: utf-8 -*-
"""
Created on Sun Apr 09 10:23:31 2017

@author: Kasia
"""

import cantera as ct
import math
import numpy as np
import matplotlib.pyplot as plt

A=0.1

#calculate kappa
def kappa(gas):
    return gas.cp / gas.cv
#return individual gas constant
def R(gas):
    return ct.gas_constant / gas.mean_molecular_weight
#calculate speed of sound
def sound(gas):
    return math.sqrt(kappa(gas) * gas.T * R(gas))
#calculate velocity from speed of sound and Mach number
def vel(M, gas):
    return M * sound(gas)
#calculate exhaust speed from inlet temperature and outlet temperature
def ve(T0, T, gas):
    return math.sqrt(2.0*gas.cp*T0*(1-T/T0))
#calculate thrust
def thrust(v, v0, mdot, gas, p0):
    return mdot*v*(v-v0)-A*(gas.P-p0)

air0 = ct.Solution('air.xml') #air
gas = ct.Solution('gri30.xml') #gas model for fuel, exhaust and igniter

#arrays for plotting
mass = []
temp = []
thr = []
thrm = []
effic = []

Mach = np.linspace(0.01,6.01,50) #Mach array

for Ma in Mach:

    #air parameters
    air = ct.Solution('air.xml')
    
    ta = 300.0
    ha = air.h
    pa = ct.one_atm
    ka = kappa(air)
    da = air.density
    
    air.TP = ta, pa
    air0.TP = ta, pa
    
    #velocity, entalphy of inlet air
    v = vel(Ma,air)
    #change entalphy because of air velocity
    h = ha + v**2.0/2.0
    #set air parameters after entalphy change
    air.HP = h, pa
    #calculate pressure
    p = pa * (air.T / ta) ** (ka/(ka-1))
    #set air parameters after pressure change
    air.HP = h, p
    
    #create reservoir for air
    air_in = ct.Reservoir(air)
    air_mw = air.mean_molecular_weight
    air_r = air.density
    
    #fuel parameters
    tf = 300.0
    pf = ct.one_atm
    xf = 'CH4:1.0' # pure methane
    
    #create reservoir for fuel
    gas.TPX = tf, pf, xf
    fuel_in = ct.Reservoir(gas)
    fuel_mw = gas.mean_molecular_weight
    fuel_r = gas.density
    
    #create reservoir for igniter
    gas.TPX = 300.0, ct.one_atm, 'H:1.0'
    igniter = ct.Reservoir(gas)
    
    #create combustion reactor
    gas.TPX = air.T, air.P, 'O2:0.21, N2:0.79'
    combustor = ct.IdealGasReactor(gas)
    combustor.volume = 1.0
    
    #create reservoir for exhaust gas
    exhaust = ct.Reservoir(gas)
    #equivalence ratio
    equiv_ratio = 1
    #stoichometric molecular air demand
    st = 9.52
    
    air_mdot = A * st * air_r
    fuel_mdot = air_mdot * equiv_ratio / st
    
    
    #create flow controllers for air and fuel
    m1 = ct.MassFlowController(fuel_in, combustor, mdot=fuel_mdot)
    m2 = ct.MassFlowController(air_in, combustor, mdot=air_mdot)
    
    #igniter time dependent mass flow
    fwhm = 0.3
    amplitude = 0.03
    t0 = 2.0
    igniter_mdot = lambda t: amplitude * math.exp(-(t-t0)**2 * 4 * math.log(2) / fwhm**2)
    m3 = ct.MassFlowController(igniter, combustor, mdot=igniter_mdot)
    
    #calculate mass flow rate
    mdot=air_mdot+fuel_mdot
    
    #create valve for combustor outlet
    valve = ct.Valve(combustor, exhaust, K=1.0)
    
    #simulation
    sim = ct.ReactorNet([combustor])
    
    t = 0.0
    tk = 4.0
    
    tempt = []
    time = []
    while (t < tk):
        t = sim.step()
        tempt.append(combustor.T)
        time.append(t)
    
    #temperature over time plotting
    """
    plt.plot(time, tempt)
    plt.xlabel('time')
    plt.ylabel('temperature')
    plt.axis([0, 4, 250, 3500])
    plt.title('Mach = %.0f' % (Ma-0.01))
    plt.savefig('M '+str('{0:.2f}'.format(Ma-0.01))+'.png')
    plt.close()
    """
    
    #calculate velocity of outlet gas
    vee = ve(combustor.T,air.T, gas)
    #calculate thrust
    th = thrust(vee,v,mdot,gas,ct.one_atm)
    #calculate engine efficiency
    eff = 2/(vee/v+1)

    mdot_v = vee * mdot
    
    #fill arrays for plotting
    if(th > 0):    
        effic.append(eff)
        temp.append(combustor.T)
        mass.append(mdot)
        thr.append(th)
        thrm.append(th/mdot_v)
    else: 
        thr.append(None)
        thrm.append(None)
        mass.append(None)
        temp.append(None)
        effic.append(None)
    
#plot results and save to .png    
names = {'Temperature': temp, 'Thrust': thr, 'Thrust over mass flow rate': thrm, 'Mass flow rate': mass, 'Efficiency': effic}

for i in names:   
    plt.plot(Mach, names[i])
    plt.xlabel('Mach Number')
    plt.ylabel('%s' % (i))
    plt.savefig(i + ' over Mach.png')
    plt.close()
