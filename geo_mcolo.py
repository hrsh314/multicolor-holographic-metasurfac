import importlib.util
# import lumapi
spec = importlib.util.spec_from_file_location('lumapi', 'C:\\Program Files\\Lumerical\\v241\\api\\python\\lumapi.py')
lumapi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lumapi)

import time
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as sc

# import custom modules
sys.path.append("../module")
from FieldPropagation import fieldPropagationLumapi, em_field
from MetaTool import nk2permittivity, setResources, integrate, phaseDis, phaseNorm

# colorbar setting
cmap_amp = "Reds"   # amplitude use
cmap_ang = "RdBu_r"   # angle (phase) use

# ============================================================================
# CONFIGURATION: WHICH SOURCE TO USE
# ============================================================================
# Set this to control which wavelength to test:
#   "red"   -> 633 nm only (test red channel)
#   "green" -> 532 nm only (test green channel)
#   "both"  -> both sources simultaneously (final multicolor test)
SOURCE_MODE = "green"  # <<< CHANGE THIS TO "green" or "both" LATER

print(f"\n{'='*70}")
print(f"  MULTICOLOR HOLOGRAM SIMULATION — SOURCE MODE: {SOURCE_MODE.upper()}")
print(f"{'='*70}\n")

# ============================================================================
# IMPORT TARGET PHASE PROFILES (20x20 each)
# ============================================================================
# Red channel phase profile (flower)
phase_profile_red = np.load("C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\phase_mask_red_20_20.npy")
# Green channel phase profile (leaves)
phase_profile_green = np.load("C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\phase_mask_green_20_20.npy")

rows, cols = phase_profile_red.shape  # should be 20x20 for both
assert phase_profile_red.shape == phase_profile_green.shape, \
    "Red and green phase profiles must have the same dimensions!"

print(f"Phase profiles loaded: {rows} × {cols}")

# Visualize target phase profiles
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

c1 = ax1.pcolor(np.rad2deg(phase_profile_red), cmap="Reds")
ax1.set_xlabel("$x$")
ax1.set_ylabel("$y$")
ax1.set_title("Target Phase — Red (633 nm, Flower)")
fig.colorbar(c1, ax=ax1)
ax1.set_aspect(1)

c2 = ax2.pcolor(np.rad2deg(phase_profile_green), cmap="Greens")
ax2.set_xlabel("$x$")
ax2.set_ylabel("$y$")
ax2.set_title("Target Phase — Green (532 nm, Leaves)")
fig.colorbar(c2, ax=ax2)
ax2.set_aspect(1)

plt.tight_layout()
plt.savefig("C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\target_phase_profiles_multicolor.png",
            dpi=300, bbox_inches='tight')
plt.show()

# ============================================================================
# PARAMETERS
# ============================================================================
# control parameters
hide = False  # whether to hide GUI or not

# spectral
wavelength_red = 633e-9    # red operating wavelength [m]
wavelength_green = 532e-9  # green operating wavelength [m]

# Choose active wavelength(s) based on SOURCE_MODE
if SOURCE_MODE == "red":
    wavelength_active = wavelength_red
    wavelength_min = wavelength_red
    wavelength_max = wavelength_red
    wavelength_number = 1
elif SOURCE_MODE == "green":
    wavelength_active = wavelength_green
    wavelength_min = wavelength_green
    wavelength_max = wavelength_green
    wavelength_number = 1
elif SOURCE_MODE == "both":
    wavelength_min = wavelength_green   # 532 nm
    wavelength_max = wavelength_red     # 633 nm
    wavelength_number = 2               # two discrete wavelengths
else:
    raise ValueError(f"Unknown SOURCE_MODE: {SOURCE_MODE}")

# ============================================================================
# UNIT CELL LAYOUT — Two meta-atoms per unit cell
# ============================================================================
#
#   Each unit cell is 400 nm × 400 nm.
#   We place TWO meta-atoms in each unit cell at diagonally opposite corners:
#
#     +-------------------+
#     |                   |
#     |  GREEN (top-left) |
#     |    ■              |
#     |                   |
#     |              ■    |
#     |   RED (bot-right) |
#     |                   |
#     +-------------------+
#
#   Corner offsets from unit cell center:
#     Red   → bottom-right corner: (+quarter, -quarter)
#     Green → top-left corner:     (-quarter, +quarter)
#
# ============================================================================

unit_size = 350e-9  # unit cell size [m] — 400 nm × 400 nm

# Offset of each meta-atom from the unit cell center
quarter = unit_size / 4  # 100 nm offset from center to each corner

units_rows = rows  # 20
units_cols = cols  # 20

# Metasurface center offsets (to center the array at origin)
offset_x = units_cols / 2 * unit_size - unit_size / 2
offset_y = units_rows / 2 * unit_size - unit_size / 2

# ============================================================================
# META-ATOM DIMENSIONS
# ============================================================================
# Red meta-atom (optimized for 633 nm)
material_metaatom_red = "Si (Silicon) - Palik"
height_red = 550e-9     # height [m] — adjust if your red library uses different
length_red = 125e-9      # x-span before rotation [m]
width_red = 95e-9       # y-span before rotation [m]

# Green meta-atom (optimized for 532 nm)
material_metaatom_green = "Si (Silicon) - Palik"
height_green = 550e-9   # height [m] — adjust if your green library uses different
length_green = 75e-9    # x-span before rotation [m]
width_green = 90e-9     # y-span before rotation [m]

# NOTE: If your red and green phase libraries were built with DIFFERENT
# meta-atom dimensions (length, width, height), update the values above
# to match each library's design parameters.

# Substrate material
material_substrate = "SiO2 (Glass) - Palik"

# Meta-atoms sit on the substrate: z = 0 to z = max(height_red, height_green)
metaatom_z_min = 0
metaatom_z_max_red = height_red
metaatom_z_max_green = height_green
metaatom_z_max = max(metaatom_z_max_red, metaatom_z_max_green)

# Substrate top at z = 0
substrate_z_max = 0

# ============================================================================
# SOURCE POSITION
# ============================================================================
source_z_position = substrate_z_max - 140e-9  # 140 nm below substrate top (inside SiO2)

# Simulation spacing — use the LONGER wavelength for safe margins
wavelength_for_spacing = wavelength_red  # 633 nm (larger of the two)
separation = wavelength_for_spacing / 2
sep_ub_t = separation
sep_t_atom = separation
sep_source_lb = separation * 0.5

# ============================================================================
# FDTD REGION Z-BOUNDARIES
# ============================================================================
fdtd_z_min = source_z_position - sep_source_lb
fdtd_z_max = metaatom_z_max + sep_t_atom + sep_ub_t

# SiO2 substrate extends past FDTD z_min into PML
substrate_z_min = fdtd_z_min - 500e-9

# Simulation size
sim_x_span = unit_size * units_cols
sim_y_span = unit_size * units_rows

# Boundary conditions
boundary_x_min = "PML"
boundary_x_max = "PML"
boundary_y_min = "PML"
boundary_y_max = "PML"
boundary_z_min = "PML"
boundary_z_max = "PML"

# Mesh settings
mesh_accuracy = 2

# ============================================================================
# Print Z-coordinate summary
# ============================================================================
print("=" * 70)
print("Z-COORDINATE LAYOUT (bottom to top)")
print("=" * 70)
print(f"  SiO2 z_min:        {substrate_z_min*1e9:.1f} nm  (extends into PML)")
print(f"  FDTD z_min:        {fdtd_z_min*1e9:.1f} nm")
print(f"  Source:             {source_z_position*1e9:.1f} nm  (140 nm below substrate top)")
print(f"  SiO2 z_max:        {substrate_z_max*1e9:.1f} nm  (substrate/meta-atom interface)")
print(f"  Meta-atom z_min:   {metaatom_z_min*1e9:.1f} nm")
print(f"  Red z_max:         {metaatom_z_max_red*1e9:.1f} nm")
print(f"  Green z_max:       {metaatom_z_max_green*1e9:.1f} nm")
print(f"  Monitor:           {(metaatom_z_max + sep_t_atom)*1e9:.1f} nm")
print(f"  FDTD z_max:        {fdtd_z_max*1e9:.1f} nm")
print("=" * 70)

# ============================================================================
# LOAD GEOMETRIC PHASE LIBRARIES
# ============================================================================

# --- RED phase library (633 nm) ---
with open('C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\633_geometric_phase_library_Si_RCP.npy', 'rb') as f:
    rotation_angle_vec_red = np.load(f)
    phase_total_vec_red = np.load(f)
    t_total_vec_red = np.load(f)

print(f"\n=== Red Phase Library (633 nm) ===")
print(f"  Rotation angles: {rotation_angle_vec_red.min():.1f}° to {rotation_angle_vec_red.max():.1f}°")
print(f"  Elements: {len(rotation_angle_vec_red)}")
print(f"  Phase range: {np.rad2deg(np.min(phase_total_vec_red)):.1f}° to {np.rad2deg(np.max(phase_total_vec_red)):.1f}°")
print(f"  Transmission range: {np.min(t_total_vec_red):.3f} to {np.max(t_total_vec_red):.3f}")

# --- GREEN phase library (532 nm) ---
with open('C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\532_geometric_phase_library_Si_RCP.npy', 'rb') as f:
    rotation_angle_vec_green = np.load(f)
    phase_total_vec_green = np.load(f)
    t_total_vec_green = np.load(f)

print(f"\n=== Green Phase Library (532 nm) ===")
print(f"  Rotation angles: {rotation_angle_vec_green.min():.1f}° to {rotation_angle_vec_green.max():.1f}°")
print(f"  Elements: {len(rotation_angle_vec_green)}")
print(f"  Phase range: {np.rad2deg(np.min(phase_total_vec_green)):.1f}° to {np.rad2deg(np.max(phase_total_vec_green)):.1f}°")
print(f"  Transmission range: {np.min(t_total_vec_green):.3f} to {np.max(t_total_vec_green):.3f}")

# Visualize both libraries
fig, axes = plt.subplots(2, 2, figsize=(14, 8))

axes[0, 0].plot(rotation_angle_vec_red, np.rad2deg(phase_total_vec_red), 'r-o', linewidth=1.5, markersize=4)
axes[0, 0].set_xlabel("Rotation Angle (deg)")
axes[0, 0].set_ylabel("Phase (deg)")
axes[0, 0].set_title("Red (633 nm) — Phase vs Rotation")
axes[0, 0].grid(True)

axes[0, 1].plot(rotation_angle_vec_red, t_total_vec_red, 'r-o', linewidth=1.5, markersize=4)
axes[0, 1].set_xlabel("Rotation Angle (deg)")
axes[0, 1].set_ylabel("Transmittance")
axes[0, 1].set_title("Red (633 nm) — Transmittance vs Rotation")
axes[0, 1].grid(True)

axes[1, 0].plot(rotation_angle_vec_green, np.rad2deg(phase_total_vec_green), 'g-o', linewidth=1.5, markersize=4)
axes[1, 0].set_xlabel("Rotation Angle (deg)")
axes[1, 0].set_ylabel("Phase (deg)")
axes[1, 0].set_title("Green (532 nm) — Phase vs Rotation")
axes[1, 0].grid(True)

axes[1, 1].plot(rotation_angle_vec_green, t_total_vec_green, 'g-o', linewidth=1.5, markersize=4)
axes[1, 1].set_xlabel("Rotation Angle (deg)")
axes[1, 1].set_ylabel("Transmittance")
axes[1, 1].set_title("Green (532 nm) — Transmittance vs Rotation")
axes[1, 1].grid(True)

plt.tight_layout()
plt.savefig("C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\phase_libraries_multicolor.png",
            dpi=300, bbox_inches='tight')
plt.show()

# ============================================================================
# FIGURE OF MERIT
# ============================================================================
def fom(phase_dest, phase, t):
    """Figure of merit for selecting best meta-atom rotation angle."""
    tolerance_phase = np.deg2rad(5)
    lambda_phase = 2e0
    lambda_t = 1e0
    return (phaseDis(phase, phase_dest) *
            (phaseDis(phase, phase_dest) > tolerance_phase) * lambda_phase +
            (-1) * t * lambda_t)

# ============================================================================
# SELECT META-ATOM ROTATION ANGLES FOR BOTH COLORS
# ============================================================================
print("\n=== Selecting Meta-atom Rotation Angles ===")

units_red_dict = {}    # Red meta-atom: position → {rotation, phase, transmission}
units_green_dict = {}  # Green meta-atom: position → {rotation, phase, transmission}

def select_rotation(phase_profile, rotation_vec, phase_vec, t_vec, color_name):
    """Select best rotation angle for each pixel from a phase library."""
    result_dict = {}
    for i in range(units_rows):
        for j in range(units_cols):
            phase_dest = phase_profile[i, j]
            matched = False
            tolerance = np.deg2rad(15)

            while not matched:
                valid_indices = np.where(np.abs(phase_vec - phase_dest) <= tolerance)[0]
                if len(valid_indices) == 0:
                    tolerance *= 2
                    if tolerance > np.pi:
                        # Fallback: pick the closest phase in the library
                        best_index = np.argmin(np.abs(phase_vec - phase_dest))
                        result_dict[(i, j)] = {
                            'rotation': rotation_vec[best_index],
                            'phase': phase_vec[best_index],
                            'transmission': t_vec[best_index]
                        }
                        matched = True
                        print(f"  [{color_name}] ({i},{j}): FALLBACK — closest phase used")
                else:
                    fom_min = np.inf
                    best_index = 0
                    for idx in valid_indices:
                        fom_current = fom(phase_dest, phase_vec[idx], t_vec[idx])
                        if fom_current < fom_min:
                            fom_min = fom_current
                            best_index = idx
                    result_dict[(i, j)] = {
                        'rotation': rotation_vec[best_index],
                        'phase': phase_vec[best_index],
                        'transmission': t_vec[best_index]
                    }
                    matched = True
                    if i < 2 and j < 2:
                        print(f"  [{color_name}] ({i},{j}): target={np.rad2deg(phase_dest):.1f}°, "
                              f"rot={rotation_vec[best_index]:.1f}°, "
                              f"phase={np.rad2deg(phase_vec[best_index]):.1f}°, "
                              f"T={t_vec[best_index]:.3f}")
    return result_dict

print("\n--- Red Channel (Flower, 633 nm) ---")
units_red_dict = select_rotation(
    phase_profile_red, rotation_angle_vec_red, phase_total_vec_red, t_total_vec_red, "RED"
)

print("\n--- Green Channel (Leaves, 532 nm) ---")
units_green_dict = select_rotation(
    phase_profile_green, rotation_angle_vec_green, phase_total_vec_green, t_total_vec_green, "GREEN"
)

# Statistics
rot_red = [units_red_dict[p]['rotation'] for p in units_red_dict]
rot_green = [units_green_dict[p]['rotation'] for p in units_green_dict]
print(f"\n=== Selection Statistics ===")
print(f"  Red:   {len(np.unique(np.round(rot_red, 1)))} unique rotations, "
      f"range [{np.min(rot_red):.1f}°, {np.max(rot_red):.1f}°]")
print(f"  Green: {len(np.unique(np.round(rot_green, 1)))} unique rotations, "
      f"range [{np.min(rot_green):.1f}°, {np.max(rot_green):.1f}°]")

# ============================================================================
# LUMERICAL FDTD SETUP
# ============================================================================
fdtd = lumapi.FDTD(hide=hide)
print("\n>> Progress: FDTD is opened.")

# Resource settings
parallel_job_number = 1
processes = 6
threads = 1
capacity = 1
job_launching_preset = "Remote: Intel MPI"

setResources(fdtd, parallel_job_number=parallel_job_number, processes=processes,
             threads=threads, capacity=capacity, job_launching_preset=job_launching_preset)

# Ensure layout mode
if fdtd.layoutmode() != 1:
    fdtd.switchtolayout()
fdtd.deleteall()

# ============================================================================
# SOURCE(S) — Controlled by SOURCE_MODE
# ============================================================================
# RCP = (Ex + i*Ey) / sqrt(2)

if SOURCE_MODE in ("red", "both"):
    # --- Red source (633 nm) ---
    fdtd.addplane(
        name="source_red_x",
        x=0, x_span=sim_x_span,
        y=0, y_span=sim_y_span,
        z=source_z_position,
        injection_axis="z", direction="forward",
        angle_theta=0, angle_phi=0,
        amplitude=1 / np.sqrt(2),
        polarization_angle=0,   # x-polarized
        phase=0,
        wavelength_start=wavelength_red,
        wavelength_stop=wavelength_red,
    )
    fdtd.addplane(
        name="source_red_y",
        x=0, x_span=sim_x_span,
        y=0, y_span=sim_y_span,
        z=source_z_position,
        injection_axis="z", direction="forward",
        angle_theta=0, angle_phi=0,
        amplitude=1 / np.sqrt(2),
        polarization_angle=90,  # y-polarized
        phase=90,               # 90° for RCP
        wavelength_start=wavelength_red,
        wavelength_stop=wavelength_red,
    )
    print(">> Progress: RED RCP source created (633 nm).")

if SOURCE_MODE in ("green", "both"):
    # --- Green source (532 nm) ---
    fdtd.addplane(
        name="source_green_x",
        x=0, x_span=sim_x_span,
        y=0, y_span=sim_y_span,
        z=source_z_position,
        injection_axis="z", direction="forward",
        angle_theta=0, angle_phi=0,
        amplitude=1 / np.sqrt(2),
        polarization_angle=0,
        phase=0,
        wavelength_start=wavelength_green,
        wavelength_stop=wavelength_green,
    )
    fdtd.addplane(
        name="source_green_y",
        x=0, x_span=sim_x_span,
        y=0, y_span=sim_y_span,
        z=source_z_position,
        injection_axis="z", direction="forward",
        angle_theta=0, angle_phi=0,
        amplitude=1 / np.sqrt(2),
        polarization_angle=90,
        phase=90,
        wavelength_start=wavelength_green,
        wavelength_stop=wavelength_green,
    )
    print(">> Progress: GREEN RCP source created (532 nm).")

# ============================================================================
# FDTD SIMULATION REGION
# ============================================================================
fdtd.addfdtd(
    dimension="3D",
    x=0.0, x_span=sim_x_span,
    y=0.0, y_span=sim_y_span,
    z_min=fdtd_z_min, z_max=fdtd_z_max,
    x_min_bc=boundary_x_min, x_max_bc=boundary_x_max,
    y_min_bc=boundary_y_min, y_max_bc=boundary_y_max,
    z_min_bc=boundary_z_min, z_max_bc=boundary_z_max,
    pml_layers=8,
    auto_shutoff_min=1e-5,
    mesh_accuracy=mesh_accuracy
)
print(">> Progress: FDTD simulation region created.")

# ============================================================================
# MONITORS
# ============================================================================
fdtd.setglobalmonitor("frequency points", wavelength_number)

# Transmission monitor (above meta-atoms)
monitor_z = metaatom_z_max + sep_t_atom
power_profile_t = fdtd.addpower(
    name="power_profile_T",
    monitor_type="2D Z-normal",
    x=0.0, x_span=sim_x_span,
    y=0.0, y_span=sim_y_span,
    z=monitor_z,
)
print(f">> Progress: Transmission monitor at z = {monitor_z*1e9:.1f} nm.")

# ============================================================================
# SUBSTRATE (SiO2) — extends past FDTD z_min into PML
# ============================================================================
fdtd.addrect(
    name="substrate_SiO2",
    x=0.0, y=0.0,
    x_span=sim_x_span, y_span=sim_y_span,
    z_min=substrate_z_min,
    z_max=substrate_z_max,
    material=material_substrate,
)
print(f">> Progress: SiO2 substrate created (z = {substrate_z_min*1e9:.1f} to {substrate_z_max*1e9:.0f} nm).")

# ============================================================================
# META-ATOMS — Two per unit cell at diagonally opposite corners
# ============================================================================
#
#   Unit cell center at (x_center, y_center).
#   Red   meta-atom → bottom-right corner: (x_center + quarter, y_center - quarter)
#   Green meta-atom → top-left corner:     (x_center - quarter, y_center + quarter)
#
print("\n=== Creating Meta-atoms (2 per unit cell) ===")
atom_count = 0

for i in range(units_rows):
    for j in range(units_cols):
        # Unit cell center
        x_center = j * unit_size - offset_x
        y_center = i * unit_size - offset_y

        # --- RED meta-atom (bottom-right corner) ---
        x_red = x_center + quarter
        y_red = y_center - quarter
        rot_red_angle = units_red_dict[(i, j)]['rotation']

        atom_name_red = f"atom_red_{i}_{j}"
        fdtd.addrect(
            name=atom_name_red,
            x=x_red, y=y_red,
            x_span=length_red,
            y_span=width_red,
            z_min=metaatom_z_min,
            z_max=metaatom_z_max_red,
            material=material_metaatom_red,
        )
        fdtd.setnamed(atom_name_red, "first axis", "z")
        fdtd.setnamed(atom_name_red, "rotation 1", rot_red_angle)
        atom_count += 1

        # --- GREEN meta-atom (top-left corner) ---
        x_green = x_center - quarter
        y_green = y_center + quarter
        rot_green_angle = units_green_dict[(i, j)]['rotation']

        atom_name_green = f"atom_green_{i}_{j}"
        fdtd.addrect(
            name=atom_name_green,
            x=x_green, y=y_green,
            x_span=length_green,
            y_span=width_green,
            z_min=metaatom_z_min,
            z_max=metaatom_z_max_green,
            material=material_metaatom_green,
        )
        fdtd.setnamed(atom_name_green, "first axis", "z")
        fdtd.setnamed(atom_name_green, "rotation 1", rot_green_angle)
        atom_count += 1

        # Print first few for verification
        if i < 2 and j < 2:
            print(f"  Unit ({i},{j}): center=({x_center*1e9:.0f},{y_center*1e9:.0f}) nm")
            print(f"    RED   @ ({x_red*1e9:.0f},{y_red*1e9:.0f}) nm, rot={rot_red_angle:.1f}°")
            print(f"    GREEN @ ({x_green*1e9:.0f},{y_green*1e9:.0f}) nm, rot={rot_green_angle:.1f}°")

print(f"\n>> Progress: Created {atom_count} meta-atoms ({atom_count//2} red + {atom_count//2} green).")

# ============================================================================
# STRUCTURE SUMMARY
# ============================================================================
print(f"\n{'='*70}")
print(f"STRUCTURE SUMMARY")
print(f"{'='*70}")
print(f"1. SiO2 substrate: z = {substrate_z_min*1e9:.1f} to {substrate_z_max*1e9:.0f} nm")
print(f"2. Source: z = {source_z_position*1e9:.0f} nm (RCP, inside SiO2)")
print(f"   Mode: {SOURCE_MODE.upper()}")
if SOURCE_MODE in ("red", "both"):
    print(f"   → Red:   633 nm")
if SOURCE_MODE in ("green", "both"):
    print(f"   → Green: 532 nm")
print(f"3. Unit cell: {unit_size*1e9:.0f} nm × {unit_size*1e9:.0f} nm, array: {units_rows}×{units_cols}")
print(f"   Red meta-atom   (bottom-right): L×W×H = {length_red*1e9:.0f}×{width_red*1e9:.0f}×{height_red*1e9:.0f} nm")
print(f"   Green meta-atom (top-left):     L×W×H = {length_green*1e9:.0f}×{width_green*1e9:.0f}×{height_green*1e9:.0f} nm")
print(f"   Corner offset from center: ±{quarter*1e9:.0f} nm")
print(f"4. Monitor: z = {monitor_z*1e9:.1f} nm")
print(f"5. FDTD: z = {fdtd_z_min*1e9:.1f} to {fdtd_z_max*1e9:.1f} nm")
print(f"6. Sim size: {sim_x_span*1e6:.1f} × {sim_y_span*1e6:.1f} μm²")
print(f"{'='*70}")

# Save simulation
save_name = f"multicolor_hologram_{SOURCE_MODE}"
save_path = f"C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\{save_name}.fsp"
fdtd.save(save_path)
print(f"\n>> Progress: Simulation saved to {save_path}")

# ============================================================================
# RUN SIMULATION
# ============================================================================
t1 = time.perf_counter()
print("\n>> Progress: Starting simulation...")
fdtd.run()
t2 = time.perf_counter()
print(f">> Progress: Simulation complete. Run time: {t2 - t1:.2f} s.")

# ============================================================================
# EXTRACT NEAR FIELD DATA
# ============================================================================
print("\n=== Extracting Near Field Data ===")

mesh_x_vec = fdtd.getdata(power_profile_t.name, 'x')[:, 0]
mesh_y_vec = fdtd.getdata(power_profile_t.name, 'y')[:, 0]

power_profile_t_e = fdtd.getresult(power_profile_t.name, 'E')
e_t_x_mat = power_profile_t_e['E'][:, :, 0, 0, 0]  # Ex
e_t_y_mat = power_profile_t_e['E'][:, :, 0, 0, 1]  # Ey
e_t_z_mat = power_profile_t_e['E'][:, :, 0, 0, 2]  # Ez

power_profile_t_h = fdtd.getresult(power_profile_t.name, 'H')
h_t_x_mat = power_profile_t_h['H'][:, :, 0, 0, 0]
h_t_y_mat = power_profile_t_h['H'][:, :, 0, 0, 1]
h_t_z_mat = power_profile_t_h['H'][:, :, 0, 0, 2]

e_t_total_magnitude = np.sqrt(np.abs(e_t_x_mat)**2 + np.abs(e_t_y_mat)**2 + np.abs(e_t_z_mat)**2)
e_t_phase = np.angle(e_t_x_mat)

t_val = fdtd.getresult(power_profile_t.name, 'T')['T'][0]

print(f"  Mesh size: {len(mesh_x_vec)} × {len(mesh_y_vec)}")
print(f"  Total transmission: {t_val:.4f}")

# ============================================================================
# VISUALIZE NEAR FIELD
# ============================================================================
sep = 1

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

im0 = axes[0].pcolormesh(mesh_x_vec[::sep]*1e6, mesh_y_vec[::sep]*1e6,
                          np.abs(e_t_x_mat[::sep, ::sep]).T, cmap=cmap_amp)
axes[0].set_xlabel(r"$x \, (\mu m)$")
axes[0].set_ylabel(r"$y \, (\mu m)$")
axes[0].set_title(rf"$|E_x|$ Amplitude — {SOURCE_MODE.upper()}")
axes[0].set_aspect('equal')
plt.colorbar(im0, ax=axes[0])

im1 = axes[1].pcolormesh(mesh_x_vec[::sep]*1e6, mesh_y_vec[::sep]*1e6,
                          np.rad2deg(np.angle(e_t_x_mat[::sep, ::sep])).T, cmap=cmap_ang)
axes[1].set_xlabel(r"$x \, (\mu m)$")
axes[1].set_ylabel(r"$y \, (\mu m)$")
axes[1].set_title(rf"$\angle E_x$ Phase (deg) — {SOURCE_MODE.upper()}")
axes[1].set_aspect('equal')
plt.colorbar(im1, ax=axes[1])

plt.tight_layout()
plt.savefig(f"C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\near_field_{SOURCE_MODE}.png",
            dpi=300, bbox_inches='tight')
plt.show()

# ============================================================================
# FAR FIELD PROPAGATION
# ============================================================================
print("\n=== Far Field Propagation ===")

observe_z = 1.0
observe_x_half_span = 1.5
observe_y_half_span = 1.5
dest_x_vec = np.linspace(-observe_x_half_span, observe_x_half_span, 161)
dest_y_vec = np.linspace(-observe_y_half_span, observe_y_half_span, 161)

print(f"  Observation plane: z = {observe_z} m")
print(f"  x: [{-observe_x_half_span}, {observe_x_half_span}] m, "
      f"y: [{-observe_y_half_span}, {observe_y_half_span}] m")
print(f"  Resolution: {len(dest_x_vec)} × {len(dest_y_vec)}")

near_field = em_field(
    [wavelength_active if SOURCE_MODE != "both" else wavelength_red],
    mesh_x_vec, mesh_y_vec,
    [monitor_z],
    fdtd.getresult(power_profile_t.name, 'E')['E'],
    fdtd.getresult(power_profile_t.name, 'H')['H']
)

fdtd.farfieldsettings("override near field mesh", True)
fdtd.farfieldsettings("near field samples per wavelength", 4)

print(">> Progress: Computing far field...")

e_far_field_observe = fieldPropagationLumapi(
    near_field, dest_x_vec, dest_y_vec,
    [observe_z], wavelength_index_vec=np.arange(0, 1), fdtd=fdtd
)

e_far_field_observe_x = e_far_field_observe[:, :, 0, 0, 0]
e_far_field_observe_y = e_far_field_observe[:, :, 0, 1, 0]
e_far_field_observe_z = e_far_field_observe[:, :, 0, 2, 0]

e_far_field_intensity_total = (np.abs(e_far_field_observe_x)**2 +
                                np.abs(e_far_field_observe_y)**2 +
                                np.abs(e_far_field_observe_z)**2)

print(">> Progress: Far field computation complete.")

# ============================================================================
# VISUALIZE FAR FIELD
# ============================================================================
# Choose colormap based on active source
if SOURCE_MODE == "red":
    far_cmap = "Reds"
    far_title = "Far Field Intensity — RED (633 nm)"
elif SOURCE_MODE == "green":
    far_cmap = "Greens"
    far_title = "Far Field Intensity — GREEN (532 nm)"
else:
    far_cmap = "hot"
    far_title = "Far Field Intensity — BOTH (633 + 532 nm)"

fig = plt.figure(figsize=(8, 7))
c = plt.pcolor(dest_x_vec, dest_y_vec, e_far_field_intensity_total.T, cmap=far_cmap)
plt.colorbar(c)
plt.xlabel(r"$x \, (m)$")
plt.ylabel(r"$y \, (m)$")
plt.title(far_title)
plt.axis("scaled")
plt.tight_layout()
plt.savefig(f"C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\far_field_{SOURCE_MODE}.png",
            dpi=300, bbox_inches='tight')
plt.show()

# ============================================================================
# EFFICIENCY CALCULATIONS
# ============================================================================
print(f"\n=== Efficiency Calculations ({SOURCE_MODE.upper()}) ===")

power_output = integrate(
    1/2 * np.sqrt(sc.epsilon_0 / sc.mu_0) *
    (np.abs(e_t_x_mat)**2 + np.abs(e_t_y_mat)**2 + np.abs(e_t_z_mat)**2),
    mesh_x_vec, mesh_y_vec
)
power_input = power_output / t_val

# Use the active wavelength for airy radius
if SOURCE_MODE == "red":
    wl_eff = wavelength_red
elif SOURCE_MODE == "green":
    wl_eff = wavelength_green
else:
    wl_eff = (wavelength_red + wavelength_green) / 2  # average for "both" mode

na = 0.6
airy_radius = 0.61 * wl_eff / na
focus_radius = 3 * airy_radius

power_focus_total = integrate(
    1/2 * np.sqrt(sc.epsilon_0 / sc.mu_0) * e_far_field_intensity_total,
    dest_x_vec, dest_y_vec, focus_radius
)

power_diffraction_total = integrate(
    1/2 * np.sqrt(sc.epsilon_0 / sc.mu_0) * e_far_field_intensity_total,
    dest_x_vec, dest_y_vec
)

eff_focal_rel = power_focus_total / power_output
eff_focal_abs = power_focus_total / power_input
eff_diff_rel = power_diffraction_total / power_output
eff_diff_abs = power_diffraction_total / power_input

print(f"\n  Transmission: {t_val:.4f} ({t_val*100:.2f}%)")
print(f"\n  Focusing Efficiency:")
print(f"    Relative: {eff_focal_rel:.4f} ({eff_focal_rel*100:.2f}%)")
print(f"    Absolute: {eff_focal_abs:.4f} ({eff_focal_abs*100:.2f}%)")
print(f"\n  Diffraction Efficiency:")
print(f"    Relative: {eff_diff_rel:.4f} ({eff_diff_rel*100:.2f}%)")
print(f"    Absolute: {eff_diff_abs:.4f} ({eff_diff_abs*100:.2f}%)")

# ============================================================================
# SAVE RESULTS
# ============================================================================
print(f"\n=== Saving Results ({SOURCE_MODE.upper()}) ===")

base_path = "C:\\Users\\HARSH\\Desktop\\Python_Lumerical"
prefix = f"multicolor_{SOURCE_MODE}"

np.save(f"{base_path}\\{prefix}_near_field_x.npy", mesh_x_vec)
np.save(f"{base_path}\\{prefix}_near_field_y.npy", mesh_y_vec)
np.save(f"{base_path}\\{prefix}_near_field_Ex.npy", e_t_x_mat)
np.save(f"{base_path}\\{prefix}_near_field_Ey.npy", e_t_y_mat)
np.save(f"{base_path}\\{prefix}_near_field_Ez.npy", e_t_z_mat)

np.save(f"{base_path}\\{prefix}_far_field_x.npy", dest_x_vec)
np.save(f"{base_path}\\{prefix}_far_field_y.npy", dest_y_vec)
np.save(f"{base_path}\\{prefix}_far_field_Ex.npy", e_far_field_observe_x)
np.save(f"{base_path}\\{prefix}_far_field_Ey.npy", e_far_field_observe_y)
np.save(f"{base_path}\\{prefix}_far_field_Ez.npy", e_far_field_observe_z)
np.save(f"{base_path}\\{prefix}_far_field_intensity.npy", e_far_field_intensity_total)

efficiency_data = {
    'source_mode': SOURCE_MODE,
    'total_transmission': t_val,
    'focal_efficiency_relative': eff_focal_rel,
    'focal_efficiency_absolute': eff_focal_abs,
    'diffraction_efficiency_relative': eff_diff_rel,
    'diffraction_efficiency_absolute': eff_diff_abs,
}
np.save(f"{base_path}\\{prefix}_efficiency_data.npy", efficiency_data)

print(">> Progress: All results saved.")
print(f"\n{'='*70}")
print(f"  SIMULATION COMPLETE — MODE: {SOURCE_MODE.upper()}")
print(f"{'='*70}")
print(f"\n  NEXT STEPS:")
if SOURCE_MODE == "red":
    print(f"  → Check if the far field shows the RED (flower) image.")
    print(f"  → If yes, change SOURCE_MODE to 'green' and re-run.")
elif SOURCE_MODE == "green":
    print(f"  → Check if the far field shows the GREEN (leaves) image.")
    print(f"  → If yes, change SOURCE_MODE to 'both' and re-run.")
elif SOURCE_MODE == "both":
    print(f"  → Check the combined multicolor far field image.")
    print(f"  → You should see both flower (red) and leaves (green)!")
