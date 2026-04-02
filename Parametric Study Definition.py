#The condtions that I can change
ratios = [4]
Tin_list = [673]
Twalls = [600]
vels = [0.1]
pressures = [7.5]

#PARAMETRIC SWEEP
results = []
profiles = []
Lz = 1.0
case = 0

#Knowing how many combinations do we have for the kinetics HOW MANY DO WE HAVE =?
total = len(ratios) * len(Twalls) * len(vels) * len(pressures) * len(Tin_list)

#The beginning of the reaction altering ratio, Tw, u and P
for r in ratios:
    for Tw in Twalls:
        for u in vels:
            for P in pressures:
                for Tin in Tin_list:
                    case += 1
                    print(f"▶ Case {case} / {total} | H2/CO2={r:.2f}, T={Tw:.2f}, u={u:.2f}, P={P:.2f}, Tin={Tin:.2f}")
    
                    Xc, Sc, Tc, Qaxialc, Tmaxc, Nr, Nz, Rr, Lz = run_reactor(r, Tw, u, P,
                                         champom_kinetics, kin_champom, "Champom", Tin)
                    Xk, Sk, Tk, Qaxialk, Tmaxk, Nr, Nz, Rr, Lz = run_reactor(r, Tw, u, P,
                                         koz_kinetics, kin_koz, "Koschany", Tin)
                    Xkai, Skai, Tkai, Qaxialkai, Tmaxkai, Nr, Nz, Rr, Lz = run_reactor(r, Tw, u, P,
                                         kai_kinetics, kin_kai, "Kai", Tin)
                    Xfar, Sfar, Tfar, Qaxialfar, Tmaxfar, Nr, Nz, Rr, Lz = run_reactor(r, Tw, u, P,
                                                     farsi_kinetics, kin_farsi, "Farsi", Tin)
                    profiles.append({ "kinetic": "Champom","H2_CO2": r,"Twall": Tw,"u0": u,"P": P,"Qaxial": Qaxialc,"Tmax": Tmaxc,"Tfield": Tc,
                                    "Nr": Nr,"Nz": Nz,"Rr": Rr,"Lz": Lz,"Tin": Tin})
                    
                    profiles.append({"kinetic": "Koz","H2_CO2": r,"Twall": Tw,"u0": u,"P": P,"Qaxial": Qaxialk,"Tmax": Tmaxk,"Tfield": Tk,
                                    "Nr": Nr,"Nz": Nz,"Rr": Rr,"Lz": Lz,"Tin": Tin})
                    
                    profiles.append({"kinetic": "Kai","H2_CO2": r,"Twall": Tw,"u0": u,"P": P,"Qaxial": Qaxialkai,"Tmax": Tmaxkai,"Tfield": Tkai,
                                    "Nr": Nr,"Nz": Nz,"Rr": Rr,"Lz": Lz,"Tin": Tin})
                    
                    profiles.append({"kinetic": "Farsi","H2_CO2": r,"Twall": Tw,"u0": u,"P": P,"Qaxial": Qaxialfar,"Tmax": Tmaxfar,"Tfield": Tfar,
                                    "Nr": Nr,"Nz": Nz,"Rr": Rr,"Lz": Lz,"Tin": Tin})
    
                    results.append([r, Tw, u, P, Xc, Sc, Xk, Sk, Xkai, Skai, Xfar, Sfar])

#Results and profiles
profiles = pd.DataFrame(profiles)
results = pd.DataFrame(
    results,
    columns=["H2_CO2", "Twall", "u0", "P",
             "X_champom", "S_champom", "X_Koz", "S_Koz","X_Kai", "S_Kai","X_far", "S_far" ]
)

results

# Add Tin to results (already in sweep loop)
results["Tin"] = Tin_list[0]  # or save per sweep if multiple


# PLOTS
# Convert conversion and selectivity to percent
results["X_ch_pct"] = results["X_champom"] * 100
results["X_koz_pct"] = results["X_Koz"] * 100
results["X_kai_pct"] = results["X_Kai"] * 100
results["X_far_pct"] = results["X_far"] * 100
results["S_ch_pct"] = results["S_champom"] * 100
results["S_koz_pct"] = results["S_Koz"] * 100
results["S_kai_pct"] = results["S_Kai"] * 100
results["S_far_pct"] = results["S_far"] * 100