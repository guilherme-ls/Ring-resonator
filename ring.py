import gdsfactory as gf
import gplugins.tidy3d as gt
import gplugins.gmeep as gm
import gplugins
import matplotlib.pyplot as plt
import numpy as np
import os
import meep as mp
from pdk import PDK, soi_strip


PDK.activate()

wavelength = 1.55
wg_width = 0.450
wg_height = 0.220
ring_radius = 5
gap_top = 0.15
gap_bot = 0.15
wg_length = 30
couple_length = 16.4787

N_SI = 3.47
N_SIO2 = 1.45

# Creates add-drop ring component
def make_add_drop_ring(radius, cross_section, gap_top, gap_bot, wg_length, couple_length):
    c = gf.Component()
    width = cross_section.width

    # Creates ring waveguide
    ring = gf.components.ring_asymmetric(radius=radius, length_x=couple_length, length_y=0, cross_section=cross_section)
    c.add_ref(ring)

    # Places straight waveguides above and below the ring
    y_top = 2 * radius + width + gap_top
    y_bot = - width - gap_bot
    x_offset = -wg_length / 2 + couple_length / 2

    # Top waveguide
    top_wg = c << gf.components.straight(length=wg_length, cross_section=cross_section)
    top_wg.move((x_offset, y_top))

    # Bottom waveguide
    bottom_wg = c << gf.components.straight(length=wg_length, cross_section=cross_section)
    bottom_wg.move((x_offset, y_bot))

    # Places ports on straight waveguides
    c.add_port(name="drop", port=top_wg.ports["o1"])
    c.add_port(name="add", port=top_wg.ports["o2"])
    c.add_port(name="input", port=bottom_wg.ports["o1"])
    c.add_port(name="through", port=bottom_wg.ports["o2"])

    return c

# Create ring component
c = make_add_drop_ring(ring_radius, soi_strip(), gap_top, gap_bot, wg_length, couple_length)

# Show layout view
c.show()
c.plot()
#scene = c.to_3d()
#scene.show()

#####################################################
# Component simulation
#####################################################

# Extract geometry from your PDK
xs = soi_strip()
layers = PDK.layer_stack.layers

core_width = xs.width
core_height = layers["core"].thickness
box_thickness = layers["box"].thickness
clad_thickness = layers["cladding"].thickness

# Straight waveguide
straight = gt.modes.Waveguide(
    wavelength=wavelength,
    core_width=core_width,
    core_thickness=core_height,
    core_material="si",
    clad_material="sio2",
    box_material="sio2",
    slab_thickness=0.0,
    box_thickness=box_thickness,
    clad_thickness=clad_thickness,
    num_modes=1,
    group_index_step=0.01,
    precision="double",
    grid_resolution=30,
)

# Bent waveguide
bent = gt.modes.Waveguide(
    wavelength=wavelength,
    core_width=core_width,
    core_thickness=core_height,
    core_material="si",
    clad_material="sio2",
    box_material="sio2",
    #core_material=3.47,
    #clad_material=1.45,
    #box_material=1.45,
    slab_thickness=0.0,
    box_thickness=box_thickness,
    clad_thickness=clad_thickness,
    bend_radius=ring_radius,
    num_modes=1,
    group_index_step=0.01,
    precision="double",
    grid_resolution=30,
)

# Extract modal properties
neff = bent.n_eff
ng = bent.n_group
loss = bent.loss_dB_per_cm
te_fraction = bent.fraction_te

neff_st = straight.n_eff
ng_st = straight.n_group

# Plot field for mode 0
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for ax, field in zip(axes, ["Ex", "Ey", "Ez"]):
    bent.plot_field(
        field_name=field,
        value="dB",
        mode_index=0,
        ax=ax,
    )
    ax.set_title(f"Mode 0: {field} (dB)")
plt.tight_layout()
plt.savefig("electric_field.png")

# Select the most TE-like mode
mode = 0
print("Effective index:", neff)
print("Group index:", ng)
print("Bend propagation loss:", loss, "dB/cm")
print("TE fraction:", te_fraction[mode], "%")
print("Effective index (straight):", neff_st)
print("Group index (straight):", ng_st)

# Calculate FSR
fsr = wavelength**2 / (2 * np.pi * ring_radius * ng + 2 * couple_length * ng_st)
print("Estimated FSR:", fsr * 1000, "nm")

####################################################
# Get effective index for 2.5D simulation
####################################################

# Simulation parameters
w_start = 1.544
w_stop = 1.556
w_point = 31
resolution = 20

# Wavelength sampling for material fitting
n_material_points = 21
wavelengths_material = np.linspace(w_start, w_stop, n_material_points)

# Silicon refractive index (Sellmeier formula)
def n_silicon(wavelength):
    l2 = wavelength**2

    n2 = (
        1
        + 10.6684293 * l2 / (l2 - 0.301516485**2)
        + 0.0030434748 * l2 / (l2 - 1.13475115**2)
        + 1.54133408 * l2 / (l2 - 1104.0**2)
    )

    return np.sqrt(n2)


# Silica refractive index (Sellmeier formula)
def n_silica(wavelength):
    l2 = wavelength**2

    n2 = (
        1
        + 0.6961663 * l2 / (l2 - 0.0684043**2)
        + 0.4079426 * l2 / (l2 - 0.1162414**2)
        + 0.8974794 * l2 / (l2 - 9.896161**2)
    )

    return np.sqrt(n2)

###########################################
# Perform a series of narrowband simulations, since a single broadband one is diverging
###########################################

wavelength_targets = np.linspace(w_start, w_stop, w_point)

half_window = 0.00025
points_per_run = 3

# Storage for the resulting S parameters
s_drop = []
s_through = []

wavelength_results = []

out_dir = "results"

os.makedirs(out_dir, exist_ok=True)

for i, wl in enumerate(wavelength_targets):

    print(f"\nSimulation {i+1}/{len(wavelength_targets)}")
    print(f"Target wavelength: {wl*1000:.3f} nm")

    # Calculate the effective core index

    n_si = n_silicon(wl)
    n_sio2 = n_silica(wl)

    n_core = gplugins.get_effective_indices(
        core_material=n_si,
        clad_materialding=n_sio2,
        nsubstrate=n_sio2,
        thickness=core_height,
        wavelength=wl,
        polarization="te",
    )[0]

    print(f"Effective core index: {n_core:.8f}")

    # Create Meep material
    core_material = mp.Medium(index=float(n_core))

    # Define a narrowband interval
    wl_min = wl - half_window
    wl_max = wl + half_window

    # Use a different output file for every wavelength
    sim_file = os.path.join(out_dir, f"ring_{wl*1000:.3f}nm.npz",)

    # Run Meep simulation
    simulation_parameters = {
        "component": c,

        "wavelength_start": wl_min,
        "wavelength_stop": wl_max,
        "wavelength_points": points_per_run,

        "resolution": resolution,
        "is_3d": False,

        "material_name_to_meep": {"si": core_material},

        "port_source_names": ["input"],

        "filepath": sim_file,
        "overwrite": True,
        "run": True,
    }

    # Parallel
    filepath = gm.write_sparameters_meep_mpi(**simulation_parameters)

    # Extract the central wavelength result
    sp = np.load(filepath)

    center_index = points_per_run // 2

    s_drop.append(sp["drop@0,input@0"][center_index])
    s_through.append(sp["through@0,input@0"][center_index])

    wavelength_results.append(wl)

    print("Simulation finished successfully.")


# Save combined results

wavelength_results = np.array(wavelength_results)
s_drop = np.array(s_drop)
s_through = np.array(s_through)

combined_file = os.path.join(out_dir, "ring_meep_nondispersive_sweep.npz")

np.savez(
    combined_file,
    wavelength=wavelength_results,
    **{
        "drop@0,input@0": s_drop,
        "through@0,input@0": s_through,
    },
)

print("\nAll simulations completed.")

# Plot

plt.figure(figsize=(9, 5))

plt.plot(wavelength_results * 1000, 20 * np.log10(np.maximum(np.abs(s_drop), 1e-15)), label="Drop")
plt.plot(wavelength_results * 1000, 20 * np.log10(np.maximum(np.abs(s_through), 1e-15)), label="Through")

plt.xlabel("Wavelength (nm)")
plt.ylabel("|S| (dB)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("s_parameters_nondispersive_sweep.png")
plt.show()