#Plot styling 
plt.rcParams.update({
    "figure.figsize": (3.35, 2.8),
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 1.5,
    "lines.markersize": 5,
    "axes.grid": False,
    "figure.dpi": 300
})

COLORS = {"Champom": "black", "Koz": "dimgray", "Kai": "silver", "Farsi": "slategray"}
MARKERS = {"Champom": "o", "Koz": ".", "Kai": ",", "Farsi": "s"}
LINESTYLES = {"Champom": "-", "Koz": "--", "Kai": ":", "Farsi": "-."}

#Generic sweep plot
def plot_sweep(ax, x, y, label):
    ax.plot(
        x, y,
        marker=MARKERS[label],
        linestyle=LINESTYLES[label],
        color=COLORS[label],
        label=label
    )

#Conversion vs Selectivity
def plot_selectivity_conversion(results):
    plt.figure(figsize=(3.35, 2.8))
    plt.scatter(results.X_ch_pct, results.S_ch_pct, marker="o", color="black", edgecolor="black", label="Champom")
    plt.scatter(results.X_koz_pct, results.S_koz_pct, marker=".",color="dimgray", edgecolor="black", label="Koz")
    plt.scatter(results.X_kai_pct, results.S_kai_pct, marker=",",color="silver", edgecolor="black", label="Kai")
    plt.scatter(results.X_far_pct, results.S_far_pct, marker="s",color="slategray", edgecolor="black", label="Farsi")
    plt.xlabel("CO$_2$ conversion / %")
    plt.ylabel("CH$_4$ selectivity / %")
    plt.title("Selectivity vs Conversion")
    plt.legend(frameon=False)
    plt.grid(True, linestyle=":", linewidth=0.5)
    plt.tight_layout()
    plt.show()

#Conversion sweeps for every parameter that was changed
def plot_conversion_sweeps(results):
    fixed_sets = results[["Twall", "u0"]].drop_duplicates().sort_values(["Twall", "u0"])
    for _, row in fixed_sets.iterrows():
        Tw, u = row.Twall, row.u0
        sub = results[(results.Twall == Tw) & (results.u0 == u)]
    
        fig, axes = plt.subplots(2, 2, figsize=(6.8, 5.2), sharey=True)
    
        # (a) Pressure
        plot_sweep(axes[0, 0], sub.P, sub.X_ch_pct, "Champom")
        plot_sweep(axes[0, 0], sub.P, sub.X_koz_pct, "Koz")
        plot_sweep(axes[0, 0], sub.P, sub.X_kai_pct, "Kai")
        plot_sweep(axes[0, 0], sub.P, sub.X_far_pct, "Farsi")
        axes[0, 0].set_xlabel("Pressure [bar]")
        axes[0, 0].set_ylabel("CO$_2$ conversion [%]")
        axes[0, 0].set_title("(a) Pressure")
    
        # (b) Wall temperature
        plot_sweep(axes[0, 1], sub.Twall, sub.X_ch_pct, "Champom")
        plot_sweep(axes[0, 1], sub.Twall, sub.X_koz_pct, "Koz")
        plot_sweep(axes[0, 1], sub.Twall, sub.X_kai_pct, "Kai")
        plot_sweep(axes[0, 1], sub.Twall, sub.X_far_pct, "Farsi")
        axes[0, 1].set_xlabel("Wall temperature [K]")
        axes[0, 1].set_title("(b) Temperature")
    
        # (c) Velocity
        plot_sweep(axes[1, 0], sub.u0, sub.X_ch_pct, "Champom")
        plot_sweep(axes[1, 0], sub.u0, sub.X_koz_pct, "Koz")
        plot_sweep(axes[1, 0], sub.u0, sub.X_kai_pct, "Kai")
        plot_sweep(axes[1, 0], sub.u0, sub.X_far_pct, "Farsi")
        axes[1, 0].set_xlabel("Superficial velocity [m s$^{-1}$]")
        axes[1, 0].set_ylabel("CO$_2$ conversion [%]")
        axes[1, 0].set_title("(c) Velocity")
    
        # (d) Ratio
        plot_sweep(axes[1, 1], sub.H2_CO2, sub.X_ch_pct, "Champom")
        plot_sweep(axes[1, 1], sub.H2_CO2, sub.X_koz_pct, "Koz")
        plot_sweep(axes[1, 1], sub.H2_CO2, sub.X_kai_pct, "Kai")
        plot_sweep(axes[1, 1], sub.H2_CO2, sub.X_far_pct, "Farsi")
        axes[1, 1].set_xlabel("H$_2$/CO$_2$ ratio [-]")
        axes[1, 1].set_title("(d) Ratio")
    
        # Shared legend
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    
        fig.suptitle(f"Twall = {Tw} K, u = {u} m s$^{{-1}}$", fontsize=9)
        plt.tight_layout(rect=[0, 0.15, 1, 1])
        plt.show()

        
#Axial heat release and Tmax for each kinetic
def plot_axial_heat_and_T(profile):

    Tin = float(profile["Tin"])
    Nr = int(profile["Nr"])
    Nz = int(profile["Nz"])
    Rr = float(profile["Rr"])
    Lz = float(profile["Lz"])

    # Safely get Qaxial or Qc
    Q = profile.get("Qaxial")
    Q = np.atleast_1d(np.array(Q, dtype=float).ravel())

    Tmax = np.atleast_1d(np.array(profile["Tmax"], dtype=float).ravel())

    # Expand scalar if needed
    if Q.size == 1 and Tmax.size > 1:
        Q = np.full_like(Tmax, Q[0])
    if Tmax.size == 1 and Q.size > 1:
        Tmax = np.full_like(Q, Tmax[0])

    if len(Q) != len(Tmax):
        # If lengths mismatch, interpolate to match
        z_Q = np.linspace(0, Lz, len(Q))
        z_T = np.linspace(0, Lz, len(Tmax))
        Q = np.interp(z_T, z_Q, Q)
        z = z_T
    else:
        z = np.linspace(0, Lz, len(Q))

    fig, ax1 = plt.subplots(figsize=(3.35, 2.8))
    ax1.plot(z, Q / 1e6, color="black", linewidth=1.5, label="Heat release")
    ax1.set_xlabel("Axial coordinate z [m]")
    ax1.set_ylabel("Volumetric heat release [MW m$^{-3}$]")

    ax2 = ax1.twinx()
    ax2.plot(z, Tmax, "--", color="firebrick", linewidth=1.5, label="T$_{max}$")
    ax2.set_ylabel("Maximum temperature [K]", color="firebrick")
    ax2.tick_params(axis="y", labelcolor="firebrick")

    # Combined legend
    lines = ax1.get_lines() + ax2.get_lines()
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, frameon=False)

    title = (
        f"{profile['kinetic']} | H2/CO2={profile['H2_CO2']}, "
        f"Twall={profile['Twall']} K, u={profile['u0']} m/s, P={profile['P']} bar"
    )
    ax1.set_title(title)
    ax1.grid(True, linestyle=":", linewidth=0.5)
    plt.tight_layout()
    plt.show()


#Radial hotspot profile
def plot_radial_hotspot(profile):
    
    Tin = float(profile["Tin"])
    Nr = int(profile["Nr"])
    Nz = int(profile["Nz"])
    Rr = float(profile["Rr"])
    Lz = float(profile["Lz"])

    Tmax = np.atleast_1d(np.array(profile["Tmax"], dtype=float).ravel())
    Tc = np.atleast_1d(np.array(profile["Tfield"], dtype=float).ravel())


    # Find axial index of hotspot
    k_hot = np.argmax(Tmax)
    Tmat = Tc.reshape(Nz, Nr)
    r = np.linspace(0, Rr, Nr)

    plt.figure(figsize=(3.35, 2.8))
    plt.plot(r * 1e3, Tmat[k_hot, :], color="black", linewidth=1.5)
    plt.xlabel("Radial coordinate r [mm]")
    plt.ylabel("Temperature [K]")
    plt.title(f"Radial profile at hotspot ({profile['kinetic']}, z ≈ {k_hot / Nz * Lz:.2f} m)")
    plt.grid(True, linestyle=":", linewidth=0.5)
    plt.tight_layout()
    plt.show()



def plot_all_profiles_combined(results, profiles):
    
    #Conversion sweeps
    plot_conversion_sweeps(results)

    #Selectivity vs conversion
    plot_selectivity_conversion(results)

    #Axial + radial profiles
    fixed_sets = profiles[["H2_CO2", "Twall", "u0", "P"]].drop_duplicates()
    
    for _, row in fixed_sets.iterrows():
        H2_CO2 = row.H2_CO2
        Twall = row.Twall
        u0 = row.u0
        P = row.P

        # Filter profiles for this parameter set
        subset = profiles[
            (profiles.H2_CO2 == H2_CO2) &
            (profiles.Twall == Twall) &
            (profiles.u0 == u0) &
            (profiles.P == P)
        ]

        # --- Combined kinetics plot ---
        plt.figure(figsize=(5, 3.5))
        ax1 = plt.gca()
        ax2 = ax1.twinx()

        for _, prof in subset.iterrows():
            kinetic = prof["kinetic"]
            Tin = float(prof["Tin"])
            Nz = int(prof["Nz"])
            Nr = int(prof["Nr"])
            Lz = float(prof["Lz"])
            Rr = float(prof["Rr"])
            
            Tmax = np.atleast_1d(np.array(prof["Tmax"]))
            Tfield = np.atleast_1d(np.array(prof["Tfield"]))
            Tavg = Tfield.reshape(Nz, Nr).mean(axis=1)
            ΔT = Tmax - Tin
            Qaxial = np.atleast_1d(np.array(prof["Qaxial"]))
            z = np.linspace(0, Lz, len(Tmax))

            ax1.plot(
                z, Tmax, '-', color=COLORS[kinetic], label=f'{kinetic} Tmax',
                marker=MARKERS[kinetic], markevery=max(Nz // 8, 1)
            )
            ax1.plot(z, Tavg, '--', color=COLORS[kinetic], label=f'{kinetic} Tavg')
            ax1.plot(z, ΔT, ':', color=COLORS[kinetic], label=f'{kinetic} ΔT')
            ax2.plot(z, Qaxial, '-.', color=COLORS[kinetic], label=f'{kinetic} Qaxial')

        ax1.set_xlabel("Axial coordinate z [m]")
        ax1.set_ylabel("Temperature [K]")
        ax2.set_ylabel("Reaction heat Q [W/m³]")
        ax1.set_title(f"H2/CO2={H2_CO2}, Twall={Twall}, u={u0}, P={P}")

        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines + lines2, labels + labels2, frameon=False, ncol=2)

        ax1.grid(True, linestyle=":", linewidth=0.5)
        plt.tight_layout()
        plt.show()

        # --- Radial hotspot profiles for each kinetic ---
        for _, prof in subset.iterrows():
            Tin = float(prof["Tin"])
            Nz = int(prof["Nz"])
            Nr = int(prof["Nr"])
            Lz = float(prof["Lz"])
            Rr = float(prof["Rr"])
            
            Tmax = np.atleast_1d(np.array(prof["Tmax"]))
            Tc = np.atleast_1d(np.array(prof["Tfield"]))

            if Tc.size != Nr * Nz:
                Tc = np.resize(Tc, (Nz * Nr))

            k_hot = np.argmax(Tmax)
            Tmat = Tc.reshape(Nz, Nr)
            r = np.linspace(0, Rr, Nr)

            plt.figure(figsize=(3.35, 2.8))
            plt.plot(r * 1e3, Tmat[k_hot, :], color=COLORS[prof["kinetic"]], linewidth=1.5,
                     label=f"{prof['kinetic']}")
            plt.xlabel("Radial coordinate r [mm]")
            plt.ylabel("Temperature [K]")
            plt.title(f"Radial profile at hotspot (z ≈ {k_hot / Nz * Lz:.2f} m)")
            plt.grid(True, linestyle=":", linewidth=0.5)
            plt.legend(frameon=False)
            plt.tight_layout()
            plt.show()



plot_all_profiles_combined(results, profiles)

def plot_two_panel_temperature_qaxial_sweep_bigger(profiles):

    # Unique combinations
    unique_sets = profiles[["H2_CO2", "Twall", "u0", "P"]].drop_duplicates()

    for _, row in unique_sets.iterrows():
        H2_CO2 = row.H2_CO2
        Twall = row.Twall
        u0 = row.u0
        P = row.P

        # Filter profiles for this parameter set
        subset = profiles[
            (profiles.H2_CO2 == H2_CO2) &
            (profiles.Twall == Twall) &
            (profiles.u0 == u0) &
            (profiles.P == P)
        ]

        # ---- Bigger figure ----
        fig, axes = plt.subplots(2, 1, figsize=(6.8, 7.0), sharex=True)

        # Secondary y-axis for Qaxial
        ax2 = axes[1].twinx()
        ax2.set_ylabel("Volumetric heat release [MW m$^{-3}$]", color="firebrick")
        ax2.tick_params(axis="y", labelcolor="firebrick")
        ax2.grid(False)

        for _, prof in subset.iterrows():
            kinetic = prof["kinetic"]
            Tmax = np.atleast_1d(np.array(prof["Tmax"]))
            Qaxial = np.atleast_1d(np.array(prof["Qaxial"]))
            Lz = float(prof["Lz"])
            Tin = float(prof["Tin"])
            Nz = len(Tmax)
            z = np.linspace(0, Lz, Nz)

            # ---- Top panel: Tmax ----
            axes[0].plot(
                z, Tmax,
                color=COLORS[kinetic],
                linestyle=LINESTYLES[kinetic],
                marker=MARKERS[kinetic],
                markevery=max(len(z)//8, 1),
                label=kinetic
            )

            # ---- Bottom panel: ΔT ----
            axes[1].plot(
                z, Tmax - Tin,
                color=COLORS[kinetic],
                linestyle=LINESTYLES[kinetic],
                marker=MARKERS[kinetic],
                markevery=max(len(z)//8, 1),
                label=f"{kinetic} ΔT"
            )

            # ---- Qaxial ----
            ax2.plot(
                z, Qaxial / 1e6,  # convert W/m³ → MW/m³
                color=COLORS[kinetic],
                linestyle='-.',
                alpha=0.7,
                label=f"{kinetic} Qaxial"
            )

        # ---------- Formatting ----------
        axes[0].set_ylabel("Maximum temperature [K]")
        axes[0].set_title(f"(a) Tmax | H2/CO2={H2_CO2}, Twall={Twall} K, u={u0} m/s, P={P} bar")

        axes[1].set_ylabel("Temperature rise $T_{max}-T_{in}$ [K]")
        axes[1].set_xlabel("Axial coordinate z [m]")
        axes[1].set_title("(b) Normalized temperature rise + Qaxial")
        axes[1].set_xlim(0, Lz)

        for ax in axes:
            ax.grid(True, linestyle=":", linewidth=0.5)

        # ---------- Legend outside ----------
        lines, labels = axes[1].get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        fig.legend(
            lines + lines2,
            labels + labels2,
            loc='upper center',
            bbox_to_anchor=(0.5, 1.02),
            ncol=2,
            frameon=False
        )

        plt.tight_layout()
        plt.show()


plot_two_panel_temperature_qaxial_sweep_bigger(profiles)


def interactive_two_panel_temperature(profiles, kinetic_filter=None):
    
    # Filter by kinetic if requested
    df = profiles.copy()
    if kinetic_filter:
        df = df[df["kinetic"] == kinetic_filter]

    # Flatten Tmax into a long DataFrame for hvplot
    rows = []
    for _, row in df.iterrows():
        Tmax_arr = np.atleast_1d(np.array(row["Tmax"]))
        Nz = len(Tmax_arr)
        z = np.linspace(0, float(row["Lz"]), Nz)
        Tin = float(row["Tin"])
        for zi, Ti in zip(z, Tmax_arr):
            rows.append({
                "z": zi,
                "Tmax": Ti,
                "DeltaT": Ti - Tin,
                "H2_CO2": row["H2_CO2"],
                "Twall": row["Twall"],
                "u0": row["u0"],
                "P": row["P"],
                "kinetic": row["kinetic"]
            })

    df_long = pd.DataFrame(rows)

    # ---------------- Panel (a): Tmax ----------------
    plot_Tmax = df_long.hvplot.line(
        x="z",
        y="Tmax",
        by="kinetic",
        width=900,
        height=400,
        hover_cols=["H2_CO2", "Twall", "u0", "P"],
        title="(a) Maximum Temperature along Reactor",
        widget_location="top"
    )

    # ---------------- Panel (b): DeltaT ----------------
    plot_dT = df_long.hvplot.line(
        x="z",
        y="DeltaT",
        by="kinetic",
        width=900,
        height=400,
        hover_cols=["H2_CO2", "Twall", "u0", "P"],
        title="(b) Temperature Rise ΔT = Tmax - Tin",
        widget_location="top"
    )

    # Stack vertically
    layout = pn.Column(plot_Tmax, plot_dT)
    return layout


interactive_two_panel_temperature(profiles)

# --- Interactive version with dropdown filters ---
def interactive_two_panel_with_filters(profiles):

    # Dropdown options
    kinetic_options = list(profiles["kinetic"].unique())
    H2_CO2_options = list(profiles["H2_CO2"].unique())
    Twall_options = list(profiles["Twall"].unique())
    u0_options = list(profiles["u0"].unique())
    P_options = list(profiles["P"].unique())

    # Panel widgets
    kinetic_w = pn.widgets.Select(name="Kinetic", options=kinetic_options, value=kinetic_options[0])
    H2_CO2_w = pn.widgets.Select(name="H2/CO2 ratio", options=H2_CO2_options, value=H2_CO2_options[0])
    Twall_w = pn.widgets.Select(name="Wall Temp [K]", options=Twall_options, value=Twall_options[0])
    u0_w = pn.widgets.Select(name="Velocity [m/s]", options=u0_options, value=u0_options[0])
    P_w = pn.widgets.Select(name="Pressure [bar]", options=P_options, value=P_options[0])

    @pn.depends(kinetic=kinetic_w, H2_CO2=H2_CO2_w, Twall=Twall_w, u0=u0_w, P=P_w)
    def update_plot(kinetic, H2_CO2, Twall, u0, P):
        # Filter the DataFrame
        df = profiles[
            (profiles["kinetic"] == kinetic) &
            (profiles["H2_CO2"] == H2_CO2) &
            (profiles["Twall"] == Twall) &
            (profiles["u0"] == u0) &
            (profiles["P"] == P)
        ]

        if df.empty:
            return pn.pane.Markdown(
                "No data available for this combination of filters.",
                style={"color": "red", "font-size": "16px"}
            )

        # Flatten Tmax
        rows = []
        for _, row in df.iterrows():
            Tmax_arr = np.atleast_1d(np.array(row["Tmax"]))
            Nz = len(Tmax_arr)
            Lz = float(row["Lz"])
            Tin = float(row["Tin"])
            z = np.linspace(0, Lz, Nz)
            for zi, Ti in zip(z, Tmax_arr):
                rows.append({
                    "z": zi,
                    "Tmax": Ti,
                    "DeltaT": Ti - Tin
                })

        df_long = pd.DataFrame(rows)

        # Plot Tmax and DeltaT
        plot_Tmax = df_long.hvplot.line(
            x="z",
            y="Tmax",
            width=900,
            height=400,
            title="Maximum Temperature along Reactor",
        )

        plot_dT = df_long.hvplot.line(
            x="z",
            y="DeltaT",
            width=900,
            height=400,
            title="Temperature Rise ΔT = Tmax - Tin",
        )

        return pn.Column(plot_Tmax, plot_dT)

    # Layout
    controls = pn.Row(kinetic_w, H2_CO2_w, Twall_w, u0_w, P_w)
    layout = pn.Column(controls, update_plot)
    return layout

# Create and serve the interactive layout
interactive_layout = interactive_two_panel_with_filters(profiles)
interactive_layout.servable()

interactive_two_panel_with_filters(profiles)