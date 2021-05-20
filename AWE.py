# -*- coding: utf-8 -*-
"""
Created on Tue May 18 10:30:10 2021

@author: user
"""

""" 
알칼라인 수전해 시스템의 Electrochemical 모델링을 위한 온도, 압력에 따른 Cell voltage를 계산하는 코드
Cell Voltage = Reversible voltage (Open-circuit voltage) + Activation overpotential + Ohmic overpotential + (Concentration overpotential)
"""

import numpy as np

class Cell_Voltage():

    def __init__(self, T, P, W, I):
               
        self.T = T
        self.P = P
        self.W = W
        self.I = I
        
        self.T0 = 298
        self.R = 8.314 # 기체상수 - 단위 (J/mol K)
        self.F = 96485.34 # Faraday 상수 - 단위 (C/mol = J/mol V)
        self.z = 2  # 반응 1몰당 전자 2개 이동 

    def CV_Abdin(self):
        
        """
        2017년 Abdin 의해 제안된 모델 
        Ref. : Z. Abdin et al., "Modelling and Simulation of an Alkaline Electrolyser Cell", Energy, 138, pp. 316-331 (2017).
        * 필요 parameters
        I: 전류[A], A: 총 전극의 표면적 [cm2], w: KOH 전해질 농도 [%], α_Ca and α_An: Electrode charge transfer coefficient,
        r_Ca and r_An: Roughness factor, i0_Ca_Ref and i0_An_Ref: reference exchange current density [A/cm2],
        G_Ca, G_An: Gibbs free energy [kJ/mol], ε_Electrode, ε_Electrolyte, ε_O2, ε_H2, ε_H2O: Porosity,
        𝛿_Ele and 𝛿_Sep: thickness [cm], l: length between electrode and seperator [cm],
        β: width of bubble zone [cm], τ: Tortuosity, ω: Wettability, S_A : Surface area of Zirfon
        """
        
        A = 300
        α_Ca = 0.73
        α_An = 1.65
        r_Ca = 1.05
        r_An = 1.25
        i0_Ca_Ref = 0.001
        i0_An_Ref = 1.00E-11
        G_Ca = 80.51
        G_An = 80.51
        ε_Electrode = 0.3
        ε_Sep = 0.42
        ε_O2 = 106.7
        ε_H2 = 59.7
        ε_H2O = 809.1
        𝛿_Ele = 0.2
        𝛿_Sep = 0.05
        l = 0.125
        β = 0.045
        τ = 3.65
        τ_Sep = 2.18
        ω = 0.85
        σ_O2 = 3.467
        σ_H2 = 2.827
        σ_H2O = 2.641
        S_A = 22000

        # i: Current density[A/cm2]
        i = self.I / A
        
        # molality: KOH molality [mol/kg]
        molality = ((self.W *1000) / 56.1056) / (1-self.W)
        # Vapour_Pressure: Vapour pressure of the KOH solution [bar]
        Vapour_Pressure = 10**(-0.01508 * molality - 0.0016788 * (molality**2) + 2.25887 * (10**-5) * (molality**3) 
                           + (1 - 0.0012062 * molality + 5.6024 * (10**-4) * (molality**2) - 7.8228 * (10**-6) * (molality**3)) 
                           * (35.4462 - (3343.93/self.T) - 10.9 * np.log10(self.T) + 0.0041645*self.T))
        # Theta: bubble coverage
        θ = (-97.25 + 182 * (self.T/self.T0) - 84 * ((self.T/self.T0)**2)) * ((i / 30)**0.3) * (self.P / (self.P-Vapour_Pressure))
        
        def Activation_Overpotential():
            # Exchange current densities [mA/cm2]
            i0_Cathode = r_Ca * np.exp((G_Ca*1000/self.R)*((1/self.T) - (1/self.T0))) * i0_Ca_Ref
            i0_Anode = r_An * np.exp((G_An*1000/self.R)*((1/self.T) - (1/self.T0))) * i0_An_Ref
            
            if self.I == 0:
                AO_An = 0
                AO_Ca = 0
                
            else:
                AO_Ca = ((self.R * self.T) / (α_Ca * self.F)) * np.log(i/(i0_Cathode*(1-θ)))
                AO_An = ((self.R * self.T) / (α_An * self.F)) * np.log(i/(i0_Anode*(1-θ)))
                
            return AO_An + AO_Ca
        
        def Ohm_Resistance():
            
            # Resistivity of the Nickel (Material) [S cm]
            ρ_Electrode = 6.84 * (10**-6)
            # Temperature coefficient of Nickel
            κ = 0.006
            
            # Resistance of Electrode [S]
            Resistance_Cathode = (ρ_Electrode / ((1-ε_Electrode)**1.5)) * (𝛿_Ele/(A)) * (1+κ*(self.T-self.T0))
            Resistance_Anode = (ρ_Electrode / ((1-ε_Electrode)**1.5)) * (𝛿_Ele/(A)) * (1+κ*(self.T-self.T0))
            
            Resistance_Electrode = Resistance_Cathode + Resistance_Anode
            
            # Ohm resistance of electrolyte [S]
            Resistance_Electrolyte = (ρ_Electrode/(1+κ*(self.T-self.T0)))*(((l-β)/A + (β/(A*((1-θ)**1.5)))*2))
            
            # Resistance of Separator [S]
            Resistance_Separator = ρ_Electrode * (τ_Sep**2) * 𝛿_Sep / (ω * ε_Sep * A)
            
            ResistanceOverpotential = (Resistance_Electrode + Resistance_Separator + Resistance_Electrolyte) * self.I
            
            return ResistanceOverpotential
        
        def Reversible_Voltage():
            # M_x = molar weight of x [kg/mol]
            M_H2O = 0.018
            M_O2 = 0.032
            M_H2= 0.002
            
            # τ = Tortuosity
            τ_O2_H2O = self.T/(ε_O2*ε_H2O)**0.5
            τ_H2_H2O = self.T/(ε_H2*ε_H2O)**0.5
            
            # DC : dimensionless diffusion collision integral
            DC_O2_H2O = 1.06/((τ_O2_H2O)**0.156) + 0.193/np.exp(0.476*τ_O2_H2O) + 1.036/np.exp(1.53*τ_O2_H2O) + 1.765/(3.894*(τ_O2_H2O))
            DC_H2_H2O = 1.06/((τ_H2_H2O)**0.156) + 0.193/np.exp(0.476*τ_H2_H2O) + 1.036/np.exp(1.53*τ_H2_H2O) + 1.765/(3.894*(τ_H2_H2O))
            
            # σ : mean molecular radii of species
            σ_O2_H2O = (σ_O2 + σ_H2O)/2
            σ_H2_H2O = (σ_H2 + σ_H2O)/2
            
            # D_eff : effective binary diffusion coefficient[m2/s]
            D_eff_O2_H2O = 0.00133*((1/M_O2 + 1/M_H2O)**0.5)*(self.T**1.5)/(self.P*(σ_O2_H2O**2)*DC_O2_H2O)
            D_eff_H2_H2O = 0.00133*((1/M_H2 + 1/M_H2O)**0.5)*(self.T**1.5)/(self.P*(σ_H2_H2O**2)*DC_H2_H2O)
            
            # ρ_B : Weight density [kg/m3]
            ρ_B = 1000
            # r : mean pore radius[m]
            r = 2*ε_Sep/(S_A*ρ_B)
                        
            D_eff_H2O = (4/3)*r*((8*self.R*self.T/(np.pi*M_H2O))**0.5)
            
            # ε_Electrode/τ : ratio of electrode porosity to tortuosity
            D_eff_an = 1/((ε_Electrode/τ)*((1/D_eff_O2_H2O) + (1/D_eff_H2O)))
            D_eff_cat = 1/((ε_Electrode/τ)*(1/D_eff_H2_H2O + 1/D_eff_H2O))
                            
            # p_x : Partial pressure of x
            p_Hy = (1/((ε_Electrode/τ)*np.exp(self.R * self.T * l * i / (2 * self.F * self.P * D_eff_cat)*(9.87*10**-4)))-1)*Vapour_Pressure
            p_Ox = (1/((ε_Electrode/τ)*np.exp(self.R * self.T * l * i / (2 * self.F * self.P * D_eff_an)*(9.87*10**-4)))-1)*Vapour_Pressure
                    
            # ΔS0 : Standard state entropy change [J/mol k]
            ΔS0 = 163.275
            
            RV = 1.23 + (self.T-self.T0) * ΔS0 / (self.z*self.F) + (self.R*self.T)/(self.z*self.F) * np.log(p_Hy*(p_Ox**0.5))
            
            return RV
            
        CV = Reversible_Voltage() + Activation_Overpotential() + Ohm_Resistance()
        
        # CV : [V]
        return CV
    

        
class Hydrogen_Production():
    
    def __init__(self, T, P, W, I):
        
        self.T = T
        self.P = P
        self.W = W
        self.I = I
        
        self.T0 = 298
        self.R = 8.314 # 기체상수 - 단위 (J/mol K)
        self.F = 96485.34 # Faraday 상수 - 단위 (C/mol = J/mol V)
        self.z = 2  # 반응 1몰당 전자 2개 이동 
                    
    def HP_Ulleberg(self):
        
        I = self.I * 1000 # A to mA
        T = self.T - 273 # K to C
        nC = 1 # number of cells 
        A = 300 # cm2
       
        if T <= 40: # Basic Model (Constant P = 5bar)
              f = 150
              f2 = 0.990
             
              η = (((I/A)**2) / (f + (I/A)**2)) * f2
              n_Hydrogen = η * (nC*self.I/(self.F*self.z))

            
        elif T > 40 and T <= 60:
            f = 200
            f2 = 0.985
            
            η = (((I/A)**2) / (f + (I/A)**2)) * f2
            n_Hydrogen = η * (nC*self.I/(self.F*self.z))
            
        else:
            f = 250
            f2 = 0.980
        
            η = (((I/A)**2) / (f + (I/A)**2)) * f2
            n_Hydrogen = η * (nC*self.I/(self.F*self.z))
        
        # n_Hydrogen : [mol/s]
        return (n_Hydrogen)*(22.4*3600)/1000
        

from AWE import Hydrogen_Production as HP   
 
class Hydrogen_Production_Efficiency(HP):
    def __init__(self, T, P, W, I):
        
        self.T = T
        self.P = P
        self.W = W
        self.I = I
        
        self.T0 = 298
        self.R = 8.314 # 기체상수 - 단위 (J/mol K)
        self.F = 96485.34 # Faraday 상수 - 단위 (C/mol = J/mol V)
        self.z = 2  # 반응 1몰당 전자 2개 이동     

    def Efficiency(self):

        Ideal_Hydrogen_Production = self.I/(2*self.F)
        Real_Hydrogen_Production = HP.HP_Ulleberg(self)
        
        Eff = Real_Hydrogen_Production / Ideal_Hydrogen_Production
        
        # Efficiency of Hydrogen production [%]
        return Eff * 100
       

