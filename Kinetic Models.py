#Kinectic Models - Champom 2019; Koschany 2016; Kai 1987; Farsi 2020
def champom_kinetics(p, T, k):
    T = numerix.clip(T, 300, 2000)
    ln10 = numerix.log(10.0)

    kco2  = k["k0co2"]  * safe_exp(-k["eaco2"] / (R*T))
    krwgs = k["k0rwgs"] * safe_exp(-k["earwgs"] / (R*T))
    kco   = k["k0co"]   * safe_exp(-k["eaco"] / (R*T))

    Keq = 137* T**(-3.998) * safe_exp(158700/(R*T))
    Kwgs = 1.0 / numerix.exp((-2.4198 + 3.855e-4 * T + 2180.9 / T) * ln10)
    Kmeth = numerix.exp((4.1002e-5 * T**2 - 0.08025 * T + 39.6039) * ln10)
    
    Kco  = k["kads0co"]  * safe_exp(-k["deltahco"]  / (R*T))
    Kh2  = k["kads0h2"]  * safe_exp(-k["deltahh2"]  / (R*T))
    Kh2o = k["kads0h2o"] * safe_exp(-k["deltahh2o"] / (R*T))
    Kco2 = k["kads0co2"] * safe_exp(-k["deltahco2"] / (R*T))

    denom = 1 + Kh2 * p["H2"] + Kco2 * p["CO2"] + Kh2o * p["H2O"] + Kco * p["CO"]

    driving1 = safe_div(1 - (p["CH4"] * p["H2O"]**2) / (p["CO2"] * p["H2"]**4 * Keq), 1)
    driving1 = numerix.clip(driving1, 0, 1)
    driving2 = safe_div(1 - (p["CO"] * p["H2O"]) / (p["CO2"] * p["H2"] * Kwgs), 1)
    driving2 = numerix.clip(driving2, 0, 1)
    driving3 = safe_div(1 - (p["CH4"] * p["H2O"]) / (p["CO"] * p["H2"]**3 * Kmeth), 1)
    driving3 = numerix.clip(driving3, 0, 1)

    r1 = kco2 * Kh2 * Kco2 * p["H2"] * p["CO2"] * driving1 / numerix.maximum(denom**2, 1e-12)
    r2 = krwgs * Kco2 * p["CO2"] * driving2 / numerix.maximum(denom, 1e-12)
    r3   = kco * Kco2 * Kh2 * p["CO"] * p["H2"] * driving3 / numerix.maximum(denom**2, 1e-12)
    
    return [r1, r2, r3]

def koz_kinetics(p, T, k):
    T = numerix.clip(T, 300, 2000)

    kf = k["k0"] * safe_exp((k["Ea"] / R) * ( (1 / 555) - (1 / T)))
    
    Keq = 137 * T**(-3.998) * safe_exp(158700 / (R * T))
    
    K_OH  = safe_exp(k["Aoh"]  + (k["Boh"] / T))
    K_H2  = safe_exp(k["Ah2"]  + (k["Bh2"] / T))
    K_mix = safe_exp(k["Amix"] + (k["Bmix"] / T))

    driving = safe_div(1 - (p["CH4"] * p["H2O"]**2)/(p["CO2"] * p["H2"]**4 * Keq),1)
    driving = numerix.clip(driving,0,1)

    num = kf * p["H2"]**0.5 * p["CO2"]**0.5 * driving
    den = (1 + K_OH*(p["H2O"] / p["H2"]**0.5) + K_H2*p["H2"]**0.5 + K_mix*p["CO2"]**0.5)**2
    r = safe_div(num, den)
    return [r, 0.0, 0.0]

def kai_kinetics(p, T, k):
    T = numerix.clip(T, 300, 2000)

    kf = k["k0_kai"] * safe_exp(-k["Ea_kai"] / (R * T))
 
    Kh2kai  = safe_exp(k["kads_h2_kai"]  + (k["deltaH_h2_kai"] / (R * T)))
    Kco2kai  = safe_exp(k["kads_co2_kai"]  + (k["deltaH_co2_kai"] / (R * T)))
    Kh2okai = safe_exp(k["kads_h2o_kai"] + (k["deltaH_h2o_kai"] / (R * T)))

    denom = 1 + Kh2kai * (p["H2"]**(1/2)) + Kco2kai * (p["CO2"]**(1/2)) + Kh2okai * p["H2O"]

    r = kf * (p["H2"]**(1/2)) * (p["CO2"]**(1/3)) / numerix.maximum(denom**2, 1e-12)

    return [r, 0.0, 0.0]

def farsi_kinetics(p, T, k):

    # Temperature limits
    T = numerix.clip(T, 300.0, 2000.0)

    # === Arrhenius constants ===
    k1 = k["k0_1"] * numerix.exp((k["ea_1"] / R) * ( (1 / 555) - (1 / T)))
    k2 = k["k0_2"] * numerix.exp((k["ea_2"] / R) * ( (1 / 555) - (1 / T)))

    # === Adsorption constant (H2O) ===
    Kh2o = k["K_H2O_Farsi"] * numerix.exp(
        (k["ea_H2O_Farsi"] / R) * ((1.0 / 555.0) - (1.0 / T))
    )

    # === Equilibrium constants (FiPy-safe) ===
    ln10 = numerix.log(10.0)

    Kwgs = 1.0 / numerix.exp(
        (-2.4198 + 3.855e-4 * T + 2180.9 / T) * ln10
    )

    Kmeth = numerix.exp(
        (4.1002e-5 * T**2 - 0.08025 * T + 39.6039) * ln10
    )

    # === Total pressure (local) ===
    Ptot = sum(p[s] for s in p)
    
    denom = 1.0 + Kh2o * p["H2O"]
    denom2 = numerix.maximum(denom**2, 1e-12)

    drivingWGS = 1.0 - (
        (p["CO"] * p["H2O"]) /
        (p["CO2"] * p["H2"]**4 * Kwgs + 1e-20)
    )
    drivingWGS = numerix.clip(drivingWGS, 0.0, 1.0)

    drivingMeth = 1.0 - (((p["CH4"] * p["H2O"]) * Ptot**2) / (p["CO"] * p["H2"]**3 * Kmeth + 1e-20))
    drivingMeth = numerix.clip(drivingMeth, 0.0, 1.0)

    r2 = k1 * (p["CO2"]**0.5) * (p["H2"]**0.5) * drivingWGS / denom2
    r3 = k2 * p["CO"] * (p["H2"]**0.5) * drivingMeth / denom2

    return [0.0, r2, r3]