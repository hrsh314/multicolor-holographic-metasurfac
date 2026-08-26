import importlib.util
spec = importlib.util.spec_from_file_location('lumapi', 'C:\\Program Files\\Lumerical\\v241\\api\\python\\lumapi.py')
lumapi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lumapi)

import time
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import klayout.db as db

# import custom modules
sys.path.append("../module")
from FieldPropagation import fieldPropagationLumapi, em_field
from MetaTool import nk2permittivity, setResources, getMatrixCenter, phaseDis

# colorbar setting
cmap_amp = "Reds"
cmap_ang = "RdBu_r"

# parameters
hide = False

# spectral
wavelength_number = 1
wavelength = 633e-9          # [m]
wavelength_min = wavelength
wavelength_max = wavelength
 
# meta-atom unit cell
unit_size = 350e-9           # [m]
offset_x = 1 / 2 * unit_size - unit_size / 2
offset_y = 1 / 2 * unit_size - unit_size / 2

# simulation objects
material_atom = "Si (Silicon) - Palik"
material_substrate = "SiO2 (Glass) - Palik"

# Layer thicknesses
height_substrate = 650e-9   # [m]
height_atom = 550e-9        # [m]

# Rectangular meta-atom dimensions
rect_width = 95e-9         # [m] Width (y-direction)
rect_length = 125e-9         # [m] Length (x-direction)

# Rotation sweep parameters
rotation_angles = np.arange(0, 181, 10)  # 0 to 180 degrees in 10-degree steps

# Source position
sep_interface_source = wavelength_max / 2 * 0.5
source_z = -sep_interface_source

separation = wavelength_max / 2
sep_ub_t = separation
sep_t_atom = separation
sep_source_lb = separation * 0.5

# simulation size
sim_x_span = unit_size
sim_y_span = unit_size
sim_z_span = height_substrate + height_atom + sep_t_atom + sep_ub_t + sep_interface_source + sep_source_lb

# boundary conditions
boundary_x_min = "Period"
boundary_x_max = "Period"
boundary_y_min = "Period"
boundary_y_max = "Period"
boundary_z_min = "PML"
boundary_z_max = "PML"

# mesh settings
mesh_accuracy = 2

# open fdtd
fdtd = lumapi.FDTD(hide=hide)
print(">> Progress: FDTD is opened.")
print(">> Progress: Using built-in materials (Si and SiO2).")

# resource settings
parallel_job_number = 6
processes = 1
threads = 1
capacity = 1
job_launching_preset = "Remote: Intel MPI"

setResources(fdtd, parallel_job_number=parallel_job_number, processes=processes,
    threads=threads, capacity=capacity, job_launching_preset=job_launching_preset)

# switch to layout mode
if fdtd.layoutmode() != 1:
    fdtd.switchtolayout()

fdtd.deleteall()

# RCP Source - X-polarized component
source_x = fdtd.addplane(
    name="source_x",
    x=0,
    x_span=sim_x_span,
    y=0,
    y_span=sim_y_span,
    z=source_z,
    injection_axis="z",
    direction="forward",
    angle_theta=0,
    angle_phi=0,
    amplitude=1 / np.sqrt(2),
    polarization_angle=0,
    phase=0,
    wavelength_start=wavelength_min,
    wavelength_stop=wavelength_max,
)

# RCP Source - Y-polarized component
source_y = fdtd.addplane(
    name="source_y",
    x=0,
    x_span=sim_x_span,
    y=0,
    y_span=sim_y_span,
    z=source_z,
    injection_axis="z",
    direction="forward",
    angle_theta=0,
    angle_phi=0,
    amplitude=1 / np.sqrt(2),
    polarization_angle=90,
    phase=90,
    wavelength_start=wavelength_min,
    wavelength_stop=wavelength_max,
)

print(">> RCP source created (two orthogonal plane sources)")

# FDTD region
sim_region = fdtd.addfdtd(
    dimension="3D",
    x=0.0,
    x_span=sim_x_span,
    y=0.0,
    y_span=sim_y_span,
    z_min=-(sep_interface_source + sep_source_lb),
    z_max=height_atom + sep_t_atom + sep_ub_t,
    x_min_bc=boundary_x_min,
    x_max_bc=boundary_x_max,
    y_min_bc=boundary_y_min,
    y_max_bc=boundary_y_max,
    z_min_bc=boundary_z_min,
    z_max_bc=boundary_z_max,
    pml_layers=8,
    auto_shutoff_min=1e-5,
    mesh_accuracy=mesh_accuracy
)

# monitor
fdtd.setglobalmonitor("frequency points", wavelength_number)
power_profile_t = fdtd.addpower(
    name="power profile T",
    monitor_type="2D Z-normal",
    x=0.0,
    x_span=sim_x_span,
    y=0.0,
    y_span=sim_y_span,
    z=height_atom + sep_t_atom,
)

# Silicon substrate
substrate = fdtd.addrect(
    name="substrate",
    x=0.0,
    y=0.0,
    x_span=sim_x_span,
    y_span=sim_y_span,
    z_max=0,
    z_min=-height_substrate,
    material=material_substrate
)

# ─────────────────────────────────────────────────────────────
# Rectangular meta-atom — placed in the BOTTOM‑RIGHT QUADRANT
# Centre at (+unit_size/4, –unit_size/4) relative to the cell centre.
# This ensures the rectangle stays completely inside the unit cell
# after rotation (half‑diagonal ≈84 nm < distance to boundaries = 100 nm).
# ─────────────────────────────────────────────────────────────
atom = fdtd.addrect(
    name="atom",
    x = +unit_size / 4,       # centre shifted right by 100 nm
    y = -unit_size / 4,       # centre shifted down by 100 nm
    x_span=rect_length,
    y_span=rect_width,
    z_min=0,
    z_max=height_atom,
    material=material_atom
)

# Results storage
phase_vec = np.zeros(len(rotation_angles))
t_vec = np.zeros(len(rotation_angles))

# Prepare and save simulation files with rotation sweep
for i in range(len(rotation_angles)):
    if fdtd.layoutmode() != 1:
        fdtd.switchtolayout()

    current_angle = rotation_angles[i]

    # Rotate meta-atom around its centre (z‑axis)
    fdtd.setnamed("atom", "first axis", "z")
    fdtd.setnamed("atom", "rotation 1", current_angle)

    file_name = "C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\geometric_phase_atom_" + str(i) + ".fsp"
    fdtd.save(file_name)
    fdtd.addjob(file_name)

fdtd.runjobs()  # run all jobs in parallel

# Obtain and process data
for i in range(len(rotation_angles)):
    file_name = "C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\geometric_phase_atom_" + str(i) + ".fsp"
    fdtd.load(file_name)

    mesh_x_vec = fdtd.getdata(power_profile_t.name, 'x').reshape(-1)
    mesh_y_vec = fdtd.getdata(power_profile_t.name, 'y').reshape(-1)

    e_x_mat = fdtd.getdata(power_profile_t.name, 'Ex')[:, :, 0, 0]
    e_y_mat = fdtd.getdata(power_profile_t.name, 'Ey')[:, :, 0, 0]

    e_net_mat = np.sqrt(np.abs(e_x_mat)**2 + np.abs(e_y_mat)**2)

    phase_vec[i] = getMatrixCenter(np.angle(e_x_mat))
    t_vec[i] = fdtd.getresult(power_profile_t.name, 'T')['T'][0]

# Plot results
fig, ax1 = plt.subplots(figsize=(7, 4))

color = 'tab:blue'
ax1.set_xlabel(r"$Rotation \enspace Angle \enspace / \enspace degrees$")
ax1.set_ylabel(r"$Transmittance$", color=color)
ax1.plot(rotation_angles, t_vec, linestyle="--", linewidth=1.5,
         marker="o", markersize=5, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel(r'$Phase \enspace shift \enspace (Net \enspace E-field)$', color=color)
ax2.plot(rotation_angles, phase_vec, linestyle="--", linewidth=1.5,
         marker="o", markersize=5, color=color)
ax2.tick_params(axis='y', labelcolor=color)


plt.title(f"Geometric Phase Library (Net E-field) — Si-Palik\n"
          f"L={rect_length*1e9:.0f} nm, W={rect_width*1e9:.0f} nm, "
          f"H={height_atom*1e9:.0f} nm, P={unit_size*1e9:.0f} nm, "
          f"λ={wavelength*1e9:.0f} nm")
plt.tight_layout()
plt.show()

# Save library
with open('C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\633_geometric_phase_final_library_Si_RCP.npy', 'wb') as f:
    np.save(f, rotation_angles)
    np.save(f, phase_vec)
    np.save(f, t_vec)

print(">> Geometric phase library saved successfully!")
print(f">> Rotation angles: {rotation_angles[0]}° to {rotation_angles[-1]}°")
print(f">> Number of samples: {len(rotation_angles)}")
