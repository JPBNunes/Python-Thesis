#Reactor
def run_reactor(H2_CO2, Twall, u0, Pbar, kinetics, kin_params, kin_name, Tin):

    #Geometry of reactor and mesh
    Rr, Lz = 0.1, 1.0
    Nr, Nz = 100, 100
    mesh = CylindricalGrid2D(dr=Rr/Nr, dz=Lz/Nz, nr=Nr, nz=Nz)
    dz = Lz / Nz
    dt = dz / u0

    #Conditions of inlet
    yCO2 = 1 / (1 + H2_CO2)
    yH2  = H2_CO2 * yCO2
    
    #Concentration at the beginning
    Cin = {
        "CO2": Pbar * yCO2 / (R * Tin),
        "H2":  Pbar * yH2  / (R * Tin),
        "CH4": 1e-8,
        "CO":  1e-8,
        "H2O": 1e-8
    }

    
    # Conc, Temp and Pressure as in the mesh
    C = {sp: CellVariable(mesh=mesh, value=Cin[sp], hasOld=True)
         for sp in species}
    
    T = CellVariable(mesh=mesh, value=Tin, hasOld=True)

    P = CellVariable(mesh=mesh, value=Pbar, hasOld=True)

    u = CellVariable(mesh=mesh, value=u0, hasOld=True)
    u_face = FaceVariable(mesh=mesh, rank=1)

    # axial velocity only → z-direction
    u_face[0] = 0.0
    u_face[1] = u.arithmeticFaceValue

    
    #BCs of Conc and Temp + Axial parameters
    for sp in species:
        C[sp].constrain(Cin[sp], mesh.facesLeft)      
        C[sp].faceGrad.constrain(0, mesh.facesRight) 
        C[sp].faceGrad.constrain(0, mesh.facesTop)   
        C[sp].faceGrad.constrain(0, mesh.facesBottom)

    T.constrain(Tin, mesh.facesLeft)
    T.constrain(Twall, mesh.facesRight)
    T.faceGrad.constrain(0, mesh.facesTop)
    T.faceGrad.constrain(0, mesh.facesBottom)

    P.constrain(Pbar, mesh.facesLeft)
    P.constrain(Pbar, mesh.facesRight)
    
    u.constrain(u0, mesh.facesLeft)
    u.constrain(u0, mesh.facesRight)

    
    Qaxial = []
    Tmax_axial = []

    #Axial starts
    for k in range(Nz):
        for sp in species:
            C[sp].updateOld()
        T.updateOld()
        P.updateOld()
        u.updateOld()
        
        #Recompute kinetics at this axial plane
        p = {s: safe(C[s] * R * T) for s in species}
        Ptot = sum(p[s] for s in species)
        rates = kinetics(p, T, kin_params)
        
      #MB of Species
        for sp in species:
            Rsp = sum(nu[sp][i] * rates[i] for i in range(3))

              #Fail-save for pressure at moment 0
            if k == 0:
                print(
                    f"[{kin_name}] "
                    f"p_H2 = {float(p['H2'].value.mean()):.4g} | "
                    f"p_CO2 = {float(p['CO2'].value.mean()):.4g} | "
                    f"p_H2O = {float(p['H2O'].value.mean()):.4g} | "
                    f"p_CO = {float(p['CO'].value.mean()):.4g} | "
                    f"p_CH4 = {float(p['CH4'].value.mean()):.4g} | "
                    f"Rate = {float(Rsp.value.mean()):.4g}"
                    )


            eq = (
                TransientTerm(coeff=eps, var=C[sp]) 
                + ConvectionTerm(coeff=u_face, var=C[sp])
                ==
                DiffusionTerm(coeff=eps * Dr, var=C[sp])
                + ImplicitSourceTerm(coeff=rho_b * Rsp, var=C[sp])
            )
            eq.solve(dt=dt)

        enforce_atomic_conservation(C, atoms)

        u_face[1] = u.arithmeticFaceValue

        update_pressure_ergun(P, C, T, u, dz)
   
        enforce_constant_pressure(C, T, Pbar, species)

        update_velocity_mass(P, C, u, Cin)

        
        #EB of Species
        Qrxn = sum(-dH[i] * rates[i] for i in range(3))

        energy = (
            TransientTerm(coeff=rhoCp, var=T)
            + ConvectionTerm(coeff=u_face * rhoCp, var=T)
            ==
            DiffusionTerm(coeff=lambda_e, var=T)
            + ImplicitSourceTerm(coeff=rho_b * Qrxn, var=T)
        )
        energy.solve(dt=dt)

        #Axial heat release (radial average at plane k)
        # reshape fields
        Qmat = Qrxn.value.reshape(Nz, Nr)
        Tmat = T.value.reshape(Nz, Nr)
        
        # axial averages
        Qaxial.append(Qmat[k, :].mean())
        Tmax_axial.append(Tmat[k, :].max())

        #Atoms check
        if k % 10 == 0:
            H = atomic_total(C, atoms, "H").value.mean()
            Cc = atomic_total(C, atoms, "C").value.mean()
            O = atomic_total(C, atoms, "O").value.mean()
        
            print(f"[z={k}] H={H:.4e}, C={Cc:.4e}, O={O:.4e}")

        #Are we done yet?
        if k % 50 == 0:
            print(f"z = {k*dz:6.3f} m | Tmax = {T.value.max():7.1f} K")
            print(
                    f"[{kin_name}] "
                    f"p_H2 = {float(p['H2'].value.min()):.4g} | "
                    f"p_CO2 = {float(p['CO2'].value.min()):.4g} | "
                    f"p_H2O = {float(p['H2O'].value.max()):.4g} | "
                    f"p_CO = {float(p['CO'].value.mean()):.4g} | "
                    f"p_CH4 = {float(p['CH4'].value.max()):.4g} | "
                    f"Rate = {float(Rsp.value.max()):.4g}"
                    )
            Ptot = sum(C[sp].value.mean() * R * T.value.mean() for sp in species)
            print(f"Total pressure = {float(Ptot):.4g}")

    #Outlet conditions and information
    CO2_out = C["CO2"].value[-Nr:].mean()
    CH4_out = C["CH4"].value[-Nr:].mean()

    #Conversion and seletivity
    X = (Cin["CO2"] - CO2_out) / Cin["CO2"]
    S = CH4_out / (Cin["CO2"] - CO2_out + 1e-12)

    return X, S, T.value, np.array(Qaxial), np.array(Tmax_axial), Nr, Nz, Rr, Lz