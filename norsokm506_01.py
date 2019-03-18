import math

def Shearstress(v_sg, v_sl , mass_g, mass_l,vol_g,vol_l,holdup,vis_g,vis_l,roughness,diameter):
    #v_sg: superficia velocity of gas m/s
    #v_sl: superficia velocity of liquid m/s
    #mass_g: mass flow of gas kg/hr
    #mass_l: mass flow of liquid kg/hr
    #vol_g: volumetric flowrate of gas m3/hr
    #vol_l: volumetric flowrate of liquid m3/hr
    #holdup: %
    #vis_g: viscosity of gas cp
    #vis_l: viscosity of liquid cp
    v_m = v_sg + v_sl
    vis_m = (vis_l * holdup + vis_g * (100 - holdup)) / 100
    #change unit from cp to n.s/m2
    vis_m = vis_m * 1000
    #calculate mixyure density, unit kg/m3
    density_g = mass_g / vol_g
    density_l = mass_l / vol_l
    density_mix = (density_l * holdup + density_g * (100 - holdup)) / 100
    friction = 0.001375 * (1 + (20000 * roughness / diameter + 10 ** 6 * vis_m / (v_m * diameter * density_mix)) ** 0.33)
    tempo = 0.5 * density_mix * friction * v_m ** 2
    return tempo


def fpH_FixT(tempe, iph):

    tempo = 7
    if tempe==5.0:
        if (iph>=3.5) and (iph<=4.6): tempo = 2.0676 + 0.2309 * iph
        if (iph>4.6) and (iph<=6.5): tempo = 4.342 - (1.061 * iph) + (0.0708 * iph ** 2)
    if tempe==15.0:
        if (iph>=3.5) and (iph<=4.6): tempo= 2.0676 - (0.2309 * iph)
        if (iph>4.6) and (iph<=6.5): tempo= 4.986 - (1.191 * iph) + (0.0708 * iph ** 2)
    if tempe==20.0:
        if (iph>=3.5) and (iph<=4.6): tempo= 2.0676 - (0.2309 * iph)           
        if (iph>4.6) and (iph<=6.5): tempo= 5.1885 - (1.2353 * iph) + (0.0708 * iph ** 2)
    if tempe==40.0:
        if (iph>=3.5) and (iph<=4.6): tempo= 2.0676 - (0.2309 * iph)           
        if (iph>4.6) and (iph<=6.5): tempo= 5.1885 - (1.2353 * iph) + (0.0708 * iph ** 2)
    if tempe==60.0:
        if (iph>=3.5) and (iph<=4.6): tempo= 1.836 - (0.1818 * iph)           
        if (iph>4.6) and (iph<=6.5): tempo= 15.444 - (6.1291 * iph) + (0.8204 * iph ** 2) - (0.0371 * iph ** 3)
    if tempe==80.0:
        if (iph>=3.5) and (iph<=4.6): tempo= 2.6727 - (0.3636 * iph)           
        if (iph>4.6) and (iph<=6.5): tempo= 331.68 * math.exp(-1.2618 * iph)
    if tempe==90.0:
        if (iph>=3.5) and (iph<=4.57): tempo= 3.1355 - (0.4673 * iph)           
        if (iph>4.57) and (iph<=5.62): tempo= 21254 * math.exp(-2.1811 * iph)
        if (iph>5.62) and (iph<=6.5): tempo= 0.4014 - (0.0538 * iph)
    if tempe==120.0:
        if (iph>=3.5) and (iph<=4.3): tempo = 1.5375 - (0.125 * iph)
        if (iph>4.3) and (iph<=5.0): tempo= 5.9757 - 1.157 * iph
        if (iph>5.0) and (iph<=6.5): tempo=  0.546125 - (0.071225 * iph)
    if tempe==150.0:
        if (iph>=3.5) and (iph<=3.8): tempo = 1
        if (iph>3.8) and (iph<=5.0): tempo= 17.634 - (7.0945 * iph) + (0.715 * iph ** 2)
        if (iph>5.0) and (iph<=6.5): tempo=  0.037
 
    return tempo

def fpH_Cal(Tempe,IpH):

    TempRange = [5.0, 15.0, 20.0, 40.0, 60.0, 80.0, 90.0, 120.0, 150.0]
    
    loc=0
    
    for i, temp_i in enumerate(TempRange):
        if temp_i > Tempe: 
            loc=i
            break
    TempLower = TempRange[loc - 1]
    TempUpper = TempRange[loc]

  
    fpHLower = fpH_FixT(TempLower, IpH)
    fpHUpper = fpH_FixT(TempUpper, IpH)

    
    tempo = (fpHUpper - fpHLower) / (TempUpper - TempLower)
    
    tempo = fpHLower + (Tempe - TempLower) * tempo
    
    return tempo
    
def pHCalculator(ITemp, IPressure, ICO2Pressure, IBicarb, IIonicStrength, CalcOfpH):

    Bicarb = IBicarb / (1000 * 61.2)
    #check if necessary
    T_Fahrenheit = ITemp * 1.8 + 32
    T_Kelvin = ITemp + 273.15
    IonicStrenght = IIonicStrength / 58.44
    P_psi = IPressure * 14.503774
    
    H0 = (14.5 / 1.00258) * 10 ** (-(2.27 + 0.00565 * T_Fahrenheit - 0.00000806 * T_Fahrenheit ** 2 + 0.075 * IonicStrenght))
    K0 = 0.00258
    k1 = 387.6 * 10 ** (-(6.41 - 0.001594 * T_Fahrenheit + 0.00000852 * T_Fahrenheit ** 2 - 0.0000307 * P_psi \
            - 0.4772 * IonicStrenght ** 0.5 + 0.118 * IonicStrenght))
    k2 = 10 ** (-(10.61 - 0.00497 * T_Fahrenheit + 0.00001331 * T_Fahrenheit ** 2 - 0.00002624 * P_psi\
        - 1.166 * IonicStrenght ** 0.5 + 0.3466 * IonicStrenght))
    Ksp_FeCO3 = 10 ** (-(10.13 + 0.0182 * ITemp))
    KW = 10 ** (-(29.3868 - 0.0737549 * T_Kelvin + 7.47881 * 10 ** (-5) * T_Kelvin ** 2))
    
    #Start calculating pH
    for pHCalcNo in range(1,CalcOfpH):
        if pHCalcNo == 1: SatpHCoeff = 0
        else: SatpHCoeff = 2 * Ksp_FeCO3 / (H0 * K0 * k1 * k2 * ICO2Pressure)
        teller = 0
        q = 0
        Hion = 10 ** (-3.5)
        if ICO2Pressure > 20: Hion = 10 ** (-2.9)
    
        while q == 0:
            fHion = SatpHCoeff * Hion ** 4 + Hion ** 3 + Bicarb * Hion ** 2 - Hion * (k1 * K0 * H0 * ICO2Pressure + KW) - 2 * k1 * k2 * K0 * H0 * ICO2Pressure
            fdHion = 4 * SatpHCoeff * Hion ** 3 + 3 * Hion * Hion + 2 * Bicarb * Hion - (k1 * K0 * H0 * ICO2Pressure + KW)
            oldHion = Hion
            Hion = oldHion - fHion / fdHion
            teller = teller + 1
            if fHion == 0: q = 1
            if teller == 100: q = 1
            if abs(Hion - oldHion) < 10 ** -6 * oldHion: q = 1
        
    
        if abs(Hion - oldHion) > 10 ** -6 * oldHion:
            if pHCalcNo == 1: 
                print("The calculation of pH did not converge")
            else:
                print("The calculation of pH in condensed water saturated with Fe2+ did not converge")
    
        if pHCalcNo == 1:
            if Hion > 0:
                tempo = round(-(math.log(Hion) / 2.3026), 2)
            else:
                tempo = 99.99
        else:
            if Hion > 0:
                tempo = round(-(math.log(Hion) / 2.3026), 2)
            else:
                tempo = 99.99
    
 
    return tempo
   
def FugacityofCO2(CO2fraction, Pressure, Temp):
        
    Temp_C = Temp + 273
    A = 10 ** (Pressure * (0.0031 - 1.4 / Temp_C))
    CO2Pressure = CO2fraction * Pressure
    tempo = A * CO2Pressure
    return tempo

def Kt(temp):
    temp_table=(5,15,20, 40, 60, 80,90, 120, 150)
    value_table=(0.42,1.59,4.762,8.927,10.695,9.949,6.250,7.770,5.203)
    for i,temp_i in enumerate(temp_table):
        if temp<temp_i:
            temp_lower=temp_table[i-1]
            temp_upper=temp_table[i]
            value_lower=value_table[i-1]
            value_upper=value_table[i]
            break
    return value_lower+ abs(value_upper-value_lower)/abs(temp_upper-temp_lower)

def Cal_Norsok(CO2fraction, Pressure, Temp, v_sg, v_sl , mass_g, mass_l,vol_g,vol_l,holdup,vis_g,vis_l,roughness,diameter, fPH_In,\
                IBicarb, IIonicStrength, CalcOfpH):
    CO2Fugacity= FugacityofCO2(CO2fraction, Pressure, Temp)
    ShearStress=Shearstress(v_sg, v_sl , mass_g, mass_l,vol_g,vol_l,holdup,vis_g,vis_l,roughness,diameter)
    ph=pHCalculator(Temp, Pressure, CO2fraction*Pressure, IBicarb, IIonicStrength, CalcOfpH)
    fPH_In=fpH_Cal(Temp, ph)
    tempo=0     
    if float(CO2fraction)!=0:
        tempo = Kt(Temp) * CO2Fugacity ** 0.62 * (ShearStress / 19) ** (0.146 + 0.0324 * math.log10(CO2Fugacity)) * fPH_In 
    return tempo


