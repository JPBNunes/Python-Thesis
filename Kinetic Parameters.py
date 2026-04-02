#Champom kinetics Parameters
kin_champom = dict(
    k0co2=1900000, eaco2=110000,
    k0rwgs=29666.66667, earwgs=97100,
    k0co=3716666.667, eaco=97300,
    kads0co=0.00239, deltahco=40600,
    kads0h2=0.000052, deltahh2=52000,
    kads0h2o=0.609, deltahh2o=14500,
    kads0co2=1.07, deltahco2=9720
)

#Koz kinetics Parameters
kin_koz = dict(
    k0=3.46e-4,
    Ea=77.5e3,
    Aoh=4.16, Boh=-2694.25,
    Ah2=-2.16, Bh2=745.73,
    Amix=-2.3, Bmix=1202.79
)

#Kai kinetics Parameters
kin_kai = dict(
    k0_kai=(9.32e3*(10**(2*(-5/6)))*10**(-3)),
    Ea_kai=72.5e3,
    kads_h2_kai=3.77e-10, deltaH_h2_kai=-90.2e3,
    kads_co2_kai=1.43e-3, deltaH_co2_kai=-29.5e3,
    kads_h2o_kai=2.75e-6, deltaH_h2o_kai=-64.3e3
)

#Farsi kinetics Parameters
kin_farsi = dict(
    k0_1=0.1425e-3, ea_1=166.55e3,
    k0_2=11.5451e-3, ea_2=60.98e3,
    K_H2O_Farsi=0.6782, ea_H2O_Farsi=11.44e3,
)