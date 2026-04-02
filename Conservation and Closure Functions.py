def update_pressure_ergun(P, C, T, u, dz):
    # Total concentration
    Ctot = sum(C[sp] for sp in C)

    # Mole fractions
    y = {sp: C[sp] / Ctot for sp in C}

    # Mean molar mass
    Mbar = sum(y[sp] * MW[sp] for sp in C)

    # Gas density
    rho = P * Mbar / (R * T)

    # Ergun pressure gradient
    dPdz = -(
        150.0 * (1 - eps)**2 / eps**3 * mu_gas * u / dp**2
        + 1.75 * (1 - eps) / eps**3 * rho * u**2 / dp
    )

    # Explicit axial march (in-place!)
    P.setValue(P.old + dPdz * dz)

    # Pressure floor (IN-PLACE)
    P.setValue(numerix.maximum(P.value, 1e3))


def atomic_total(C, atoms, element):
    total = 0.0
    for sp in C:
        n = atoms[sp][element]
        if n > 0:
            total += n * C[sp]
    return total

def enforce_atomic_conservation(C, atoms, eps=1e-20):

    # Reference totals (before correction)
    H0 = atomic_total(C, atoms, "H")
    C0 = atomic_total(C, atoms, "C")
    O0 = atomic_total(C, atoms, "O")

    # Current totals
    H1 = atomic_total(C, atoms, "H")
    C1 = atomic_total(C, atoms, "C")
    O1 = atomic_total(C, atoms, "O")

    # Correction factors (cellwise)
    fH = H0 / (H1 + eps)
    fC = C0 / (C1 + eps)
    fO = O0 / (O1 + eps)

    # Apply corrections
    for sp in C:
        aH = atoms[sp]["H"]
        aC = atoms[sp]["C"]
        aO = atoms[sp]["O"]

        correction = (
            (fH ** (aH / 2.0)) *
            (fC ** aC) *
            (fO ** aO)
        )

        C[sp].value[:] *= correction.value

        # Positivity guard
        C[sp].value[:] = numerix.maximum(C[sp].value, eps)


def enforce_constant_pressure(C, T, Pbar, species):
    Ctot = sum(C[sp] for sp in species)
    Ctot_ref = Pbar / (R * T)
    factor = Ctot_ref / (Ctot + 1e-30)

    for sp in species:
        C[sp].setValue(C[sp] * factor)


def update_velocity_mass(P, C, u, C_inlet):

    # Total mass concentration (rho = sum_i C_i * M_i)
    rho = sum(C[sp].value * MW[sp] for sp in species)      # kg/m^3 at current z
    rho_in = sum(C_inlet[sp] * MW[sp] for sp in species)  # kg/m^3 at inlet

    # Update velocity to conserve mass
    u.value[:] = u.value[0] * (rho_in / rho)