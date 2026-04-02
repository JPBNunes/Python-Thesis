#Species and their values
species = ["H2", "CO2", "CH4", "H2O", "CO"]

nu = {
    "H2":  [-4, -1, -3],
    "CO2": [-1, -1,  0],
    "CH4": [ 1,  0,  1],
    "H2O": [ 2,  1,  1],
    "CO":  [ 0,  1, -1]
}

#Atomic composition per species (atoms per molecule)
atoms = {
    "H2":  {"H": 2, "C": 0, "O": 0},
    "CO2": {"H": 0, "C": 1, "O": 2},
    "CO":  {"H": 0, "C": 1, "O": 1},
    "CH4": {"H": 4, "C": 1, "O": 0},
    "H2O": {"H": 2, "C": 0, "O": 1},
}

#Constants / Parameters
R = 8.314
rho_b = 2450
eps = 0.4
Dr = 1e-5
lambda_e = 0.8
dH = [-165e3, 41e3, -206e3]
rhoCp = 1.5e6            
dp  = 4.0e-3           
A   = 0.01              
mu_gas = 2.5e-5  #Viscosity [Pa·s]

#Molar masses [kg/mol]
MW = {
    "H2": 2.016e-3,
    "CO2": 44.01e-3,
    "CO": 28.01e-3,
    "CH4": 16.04e-3,
    "H2O": 18.015e-3
}