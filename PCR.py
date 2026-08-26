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

# Try to import simpson (newer) or simps (older) for integration
try:
    from scipy.integrate import simpson
    # If simpson is available, use it (newer scipy versions)
    simps = simpson
except ImportError:
    try:
        from scipy.integrate import simps
        # If simps is available, use it (older scipy versions)
    except ImportError:
        # If neither is available, define a simple trapezoidal integration
        def simps(y, x=None, dx=1.0, axis=-1):
            if x is None:
                return np.trapz(y, dx=dx, axis=axis)
            else:
                return np.trapz(y, x=x, axis=axis)

import os

# import custom modules
sys.path.append("../module")
from FieldPropagation import fieldPropagationLumapi, em_field
from MetaTool import nk2permittivity, setResources, getMatrixCenter, phaseDis

# colorbar setting
cmap_amp = "Reds"  # amplitude use
cmap_ang = "RdBu_r"  # angle (phase) use

# parameters
# control parameters
hide = False  # whether to hide GUI or not

# spectral
wavelength_number = 1  # the number of discrete points of the spectral
wavelength = 633e-9  # [m]
wavelength_min = wavelength
wavelength_max = wavelength
source_polarization = np.deg2rad(0)  # the angle of polarization to the x axis [rad]

# meta-atom
unit_size = 400e-9  # unit size [m]

# Rectangular meta-atom parameters - NESTED PARAMETER SWEEP
width_vec = np.arange(50e-9, 200e-9, 10e-9)  # width (x-direction) [m]
length_vec = np.arange(50e-9, 200e-9, 10e-9)  # length (y-direction) [m]

# Create meshgrid for 2D parameter sweep
width_mesh, length_mesh = np.meshgrid(width_vec, length_vec, indexing='ij')

# simulation objects
material_atom = "Si (Silicon) - Palik"
material_substrate = "SiO2 (Glass) - Palik"

# Layer thicknesses
height_SiO2 = 650e-9  # [m]
height_atom = 350e-9  # [m] - Si nano-pillar height

# Source position (inside SiO2, 130 nm below meta-atom)
source_z = -130e-9  # relative to substrate/atom interface at z=0

separation = wavelength_max / 2  # safe spacing between the objects and simulation boundaries
sep_ub_t = separation  # spacing between upper bound and transmission plane
sep_t_atom = separation  # spacing between transmission plane and atom
sep_interface_source = separation * 0.5  # spacing between interface (atom / substrate) and source
sep_source_lb = separation * 0.5  # spacing between source and lower bound

# simulation size
sim_x_span = unit_size
sim_y_span = unit_size
sim_z_span = height_SiO2 + height_atom + sep_t_atom + sep_ub_t + sep_source_lb

# boundary conditions: PML / Period / Bloch / (Anti-)Symmetric
boundary_x_min = "Period"
boundary_x_max = "Period"
boundary_y_min = "Period"
boundary_y_max = "Period"
boundary_z_min = "PML"
boundary_z_max = "PML"

# mesh settings (automate mesh)
mesh_accuracy = 2

# open fdtd
fdtd = lumapi.FDTD(hide=hide)
print(">> Progress: FDTD is opened.")

# resource settings 
parallel_job_number = 6
processes = 1
threads = 1
capacity = 1
job_launching_preset = "Remote: Intel MPI"  # "Remote: Microsoft MPI" / "Remote: Intel MPI"

setResources(fdtd, parallel_job_number=parallel_job_number, processes=processes, \
    threads=threads, capacity=capacity, job_launching_preset=job_launching_preset)

# Initialize 2D arrays for results
phase_mat = np.zeros((len(width_vec), len(length_vec)))
diff_eff_mat = np.zeros((len(width_vec), len(length_vec)))  # Diffraction efficiency matrix

# Nested parameter sweep
print(f">> Progress: Starting nested parameter sweep")
print(f">> Width range: {width_vec[0]*1e9:.0f} nm to {width_vec[-1]*1e9:.0f} nm, {len(width_vec)} points")
print(f">> Length range: {length_vec[0]*1e9:.0f} nm to {length_vec[-1]*1e9:.0f} nm, {len(length_vec)} points")
print(f">> Total simulations: {len(width_vec) * len(length_vec)}")

# Flatten the parameter grid for parallel processing
simulation_count = 0
for i, width in enumerate(width_vec):
    for j, length in enumerate(length_vec):
        # switch layout
        if fdtd.layoutmode() != 1:
            fdtd.switchtolayout()
        
        fdtd.deleteall()  # clear objects
        
        # Source - position inside SiO2 substrate
        source = fdtd.addplane(
            name="source",
            # size
            x=0,
            x_span=sim_x_span,
            y=0,
            y_span=sim_y_span,
            z=source_z,  # -130 nm relative to substrate/atom interface
            # propagation direction
            injection_axis="z",
            direction="forward",
            angle_theta=0,
            angle_phi=0,
            amplitude=1,
            # polarization direction
            polarization_angle=np.rad2deg(source_polarization),
            # phase
            phase=0,
            # bandwidth
            wavelength_start=wavelength_min,
            wavelength_stop=wavelength_max,
        )

        # FDTD simulation region
        sim_region = fdtd.addfdtd(
            dimension="3D",
            x=0.0,
            x_span=sim_x_span,
            y=0.0,
            y_span=sim_y_span,
            z_min=-(sep_interface_source + sep_source_lb),
            z_max=height_atom + sep_t_atom + sep_ub_t,
            # boundary condition
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

        # Monitor for near-field (transmission plane)
        fdtd.setglobalmonitor("frequency points", wavelength_number)  # global settings
        power_profile_t = fdtd.addpower(
            name="power profile T",
            monitor_type="2D Z-normal",
            x=0.0,
            x_span=sim_x_span,
            y=0.0,
            y_span=sim_y_span,
            z=height_atom + sep_t_atom,  # above Si nano-pillar
        )
        
        # Structure
        # SiO2 substrate
        substrate = fdtd.addrect(
            name="substrate",
            x=0.0, 
            y=0.0,
            x_span=sim_x_span,
            y_span=sim_y_span,
            z_max=0,  # substrate/atom interface at z=0
            z_min=-height_SiO2,
            material=material_substrate
        )

        # Si meta-atom - Rectangular nano-pillar
        atom = fdtd.addrect(
            name="atom",
            x=0.0,
            y=0.0,
            x_span=width,  # width in x-direction
            y_span=length,  # length in y-direction
            z_min=0,  # On top of substrate
            z_max=height_atom,  # Si nano-pillar height
            material=material_atom
        )

        # save simulation with unique index
        file_name = f"C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\rect_atom_{i}_{j}.fsp"
        fdtd.save(file_name)
        fdtd.addjob(file_name)
        
        simulation_count += 1
        if simulation_count % 10 == 0:
            print(f"  Prepared {simulation_count}/{len(width_vec)*len(length_vec)} simulations")

print(">> Progress: Running jobs in parallel...")
fdtd.runjobs()  # run jobs in parallel

print(">> Progress: Collecting simulation results...")
for i, width in enumerate(width_vec):
    for j, length in enumerate(length_vec):
        file_name = f"C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\rect_atom_{i}_{j}.fsp"
        
        try:
            fdtd.load(file_name)
            
            # Get near-field data
            mesh_x_vec = fdtd.getdata(power_profile_t.name, 'x').flatten()
            mesh_y_vec = fdtd.getdata(power_profile_t.name, 'y').flatten()
            e_x_mat = fdtd.getdata(power_profile_t.name, 'Ex')[:, :, 0, 0]
            phase_mat[i, j] = getMatrixCenter(np.angle(e_x_mat))
            
            # Get full E-field and H-field data from the transmission monitor
            E_data = fdtd.getresult(power_profile_t.name, 'E')
            H_data = fdtd.getresult(power_profile_t.name, 'H')
            
            # Extract field components
            e_t_x_mat = E_data['E'][:, :, 0, 0, 0]
            e_t_y_mat = E_data['E'][:, :, 0, 0, 1]
            e_t_z_mat = E_data['E'][:, :, 0, 0, 2]
            
            h_t_x_mat = H_data['H'][:, :, 0, 0, 0]
            h_t_y_mat = H_data['H'][:, :, 0, 0, 1]
            h_t_z_mat = H_data['H'][:, :, 0, 0, 2]
            
            # Calculate Poynting vector (power density) in z-direction
            # P = 0.5 * Re(E × H*)
            # For z-component (propagation direction):
            P_z = 0.5 * np.real(e_t_x_mat * np.conj(h_t_y_mat) - e_t_y_mat * np.conj(h_t_x_mat))
            
            # Calculate total transmitted power by integrating P_z over the monitor
            if len(mesh_x_vec) > 1 and len(mesh_y_vec) > 1:
                dx = mesh_x_vec[1] - mesh_x_vec[0]
                dy = mesh_y_vec[1] - mesh_y_vec[0]
                
                # Use 2D integration (trapezoidal rule if simps not available)
                try:
                    # Try Simpson's rule
                    power_output = simps(simps(P_z, mesh_x_vec, axis=0), mesh_y_vec)
                except:
                    # Fall back to trapezoidal rule
                    power_output = np.trapz(np.trapz(P_z, mesh_x_vec, axis=0), mesh_y_vec)
            else:
                # Simple integration if mesh vectors are not properly defined
                dx = sim_x_span / len(mesh_x_vec) if len(mesh_x_vec) > 0 else sim_x_span
                dy = sim_y_span / len(mesh_y_vec) if len(mesh_y_vec) > 0 else sim_y_span
                power_output = np.sum(P_z) * dx * dy
            
            # Calculate incident power (source amplitude = 1)
            # For a plane wave with E0 = 1 V/m in free space:
            incident_intensity = 0.5 * np.sqrt(sc.epsilon_0 / sc.mu_0)  # |E| = 1
            incident_power = incident_intensity * (sim_x_span * sim_y_span)
            
            # Total transmittance (for reference)
            total_transmittance = power_output / incident_power
            
            # =============================================
            # CALCULATE DIFFRACTION EFFICIENCY
            # =============================================
            
            # Method: Use Fourier transform of near-field
            # Calculate the spatial frequency spectrum from total E-field
            E_total = np.sqrt(np.abs(e_t_x_mat)**2 + np.abs(e_t_y_mat)**2 + np.abs(e_t_z_mat)**2)
            
            # 2D Fourier transform to get diffraction pattern
            E_fft = np.fft.fftshift(np.fft.fft2(E_total))
            E_fft_power = np.abs(E_fft)**2
            
            # Normalize to incident power
            E_fft_power_normalized = E_fft_power / (incident_power * len(mesh_x_vec) * len(mesh_y_vec))
            
            # Find the 0th order (DC component) - at the center of the FFT
            center_x = len(mesh_x_vec) // 2
            center_y = len(mesh_y_vec) // 2
            
            # Power in 0th order (specular transmission)
            power_0th = E_fft_power_normalized[center_y, center_x] if center_x < len(mesh_x_vec) and center_y < len(mesh_y_vec) else 0
            
            # Method 2: Direct calculation from average Poynting vector
            avg_power_density = np.mean(P_z)
            power_0th_direct = avg_power_density * (sim_x_span * sim_y_span) / incident_power
            
            # Use the better estimate
            diff_eff_0th = max(power_0th, power_0th_direct)
            
            # Store in matrix
            diff_eff_mat[i, j] = diff_eff_0th
            
            if (i * len(length_vec) + j) % 10 == 0:
                print(f"  Processed {i*len(length_vec) + j + 1}/{len(width_vec)*len(length_vec)}: "
                      f"Width={width*1e9:.0f}nm, Length={length*1e9:.0f}nm, "
                      f"Diff Eff={diff_eff_0th*100:.1f}%, Phase={phase_mat[i, j]:.3f} rad")
                      
        except Exception as e:
            print(f"  Error processing {file_name}: {e}")
            phase_mat[i, j] = np.nan
            diff_eff_mat[i, j] = np.nan

# =============================================
# PLOT 1: Heatmap of Phase Shift
# =============================================
fig1, ax1 = plt.subplots(figsize=(10, 8))

im1 = ax1.imshow(phase_mat, 
                 extent=[length_vec[0]*1e9, length_vec[-1]*1e9, 
                         width_vec[0]*1e9, width_vec[-1]*1e9],
                 aspect='auto', 
                 cmap='rainbow',
                 origin='lower',
                 interpolation='nearest')
ax1.set_xlabel('Length (y-direction, nm)', fontsize=12)
ax1.set_ylabel('Width (x-direction, nm)', fontsize=12)
ax1.set_title('Phase Shift of Rectangular Si Meta-Atom (rad)', fontsize=14, fontweight='bold')
ax1.tick_params(axis='both', which='major', labelsize=10)
cbar1 = fig1.colorbar(im1, ax=ax1)
cbar1.set_label('Phase (rad)', fontsize=12)
cbar1.ax.tick_params(labelsize=10)

# Save Plot 1
plot1_path = 'C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\phase_heatmap.png'
fig1.savefig(plot1_path, dpi=300, bbox_inches='tight')
print(f">> Saved Plot 1: {plot1_path}")
plt.close(fig1)  # Close the figure to free memory

# =============================================
# PLOT 2: Heatmap of Diffraction Efficiency
# =============================================
fig2, ax2 = plt.subplots(figsize=(10, 8))

# Filter out NaN values for plotting
diff_eff_mat_plot = np.copy(diff_eff_mat)
diff_eff_mat_plot[np.isnan(diff_eff_mat_plot)] = 0

im2 = ax2.imshow(diff_eff_mat_plot * 100,  # Convert to percentage
                 extent=[length_vec[0]*1e9, length_vec[-1]*1e9, 
                         width_vec[0]*1e9, width_vec[-1]*1e9],
                 aspect='auto', 
                 cmap='viridis',
                 origin='lower',
                 interpolation='nearest')

ax2.set_xlabel('Length (y-direction, nm)', fontsize=12)
ax2.set_ylabel('Width (x-direction, nm)', fontsize=12)
ax2.set_title('0th Order Diffraction Efficiency of Rectangular Si Meta-Atom (%)', fontsize=14, fontweight='bold')
ax2.tick_params(axis='both', which='major', labelsize=10)
cbar2 = fig2.colorbar(im2, ax=ax2)
cbar2.set_label('Efficiency (%)', fontsize=12)
cbar2.ax.tick_params(labelsize=10)

# Save Plot 2
plot2_path = 'C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\diffraction_efficiency_heatmap.png'
fig2.savefig(plot2_path, dpi=300, bbox_inches='tight')
print(f">> Saved Plot 2: {plot2_path}")
plt.close(fig2)  # Close the figure to free memory

# =============================================
# PLOT 3: Scatter plot of Phase vs Aspect Ratio colored by Diffraction Efficiency
# =============================================
aspect_ratio_mat = width_mesh / length_mesh

fig3, ax3 = plt.subplots(figsize=(10, 8))

# Flatten arrays for scatter plot
aspect_flat = aspect_ratio_mat.flatten()
phase_flat = phase_mat.flatten()
diff_eff_flat = diff_eff_mat.flatten()

# Remove NaN values
valid_idx = ~np.isnan(phase_flat) & ~np.isnan(diff_eff_flat) & ~np.isnan(aspect_flat)
aspect_flat_valid = aspect_flat[valid_idx]
phase_flat_valid = phase_flat[valid_idx]
diff_eff_flat_valid = diff_eff_flat[valid_idx]

if len(aspect_flat_valid) > 0:
    sc3 = ax3.scatter(aspect_flat_valid, phase_flat_valid, c=diff_eff_flat_valid * 100, 
                     cmap='viridis', s=50, alpha=0.7, edgecolors='black', 
                     linewidth=0.5)

    ax3.set_xlabel('Aspect Ratio (Width/Length)', fontsize=12)
    ax3.set_ylabel('Phase Shift (rad)', fontsize=12)
    ax3.set_title('Phase vs Aspect Ratio (colored by 0th Order Diffraction Efficiency)', 
                  fontsize=14, fontweight='bold')
    ax3.tick_params(axis='both', which='major', labelsize=10)
    ax3.grid(True, alpha=0.3)
    cbar3 = fig3.colorbar(sc3, ax=ax3)
    cbar3.set_label('Diffraction Efficiency (%)', fontsize=12)
    cbar3.ax.tick_params(labelsize=10)
else:
    ax3.text(0.5, 0.5, 'No valid data to plot', 
             horizontalalignment='center', verticalalignment='center',
             transform=ax3.transAxes, fontsize=14)
    ax3.set_title('Phase vs Aspect Ratio', fontsize=14, fontweight='bold')

# Save Plot 3
plot3_path = 'C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\phase_vs_aspect_ratio.png'
fig3.savefig(plot3_path, dpi=300, bbox_inches='tight')
print(f">> Saved Plot 3: {plot3_path}")
plt.close(fig3)  # Close the figure to free memory

# =============================================
# PLOT 4: Scatter plot of Diffraction Efficiency vs Aspect Ratio
# =============================================
fig4, ax4 = plt.subplots(figsize=(10, 8))

if len(aspect_flat_valid) > 0:
    # Calculate mean diffraction efficiency for each aspect ratio
    unique_aspects = np.unique(np.round(aspect_flat_valid, 2))
    mean_efficiencies = []
    std_efficiencies = []
    for aspect in unique_aspects:
        mask = np.abs(aspect_flat_valid - aspect) < 0.01
        if np.sum(mask) > 0:
            mean_efficiencies.append(np.mean(diff_eff_flat_valid[mask] * 100))
            std_efficiencies.append(np.std(diff_eff_flat_valid[mask] * 100))

    if len(unique_aspects) > 0 and len(mean_efficiencies) > 0:
        # Trim to match lengths
        min_len = min(len(unique_aspects), len(mean_efficiencies), len(std_efficiencies))
        unique_aspects = unique_aspects[:min_len]
        mean_efficiencies = mean_efficiencies[:min_len]
        std_efficiencies = std_efficiencies[:min_len]
        
        # Plot individual points with some transparency
        ax4.scatter(aspect_flat_valid, diff_eff_flat_valid * 100, 
                   alpha=0.3, s=20, color='blue', label='Individual simulations')

        # Sort by aspect ratio
        sort_idx = np.argsort(unique_aspects)
        unique_aspects_sorted = unique_aspects[sort_idx]
        mean_eff_sorted = np.array(mean_efficiencies)[sort_idx]
        std_eff_sorted = np.array(std_efficiencies)[sort_idx]
        
        # Plot mean line
        ax4.plot(unique_aspects_sorted, mean_eff_sorted, 'r-', linewidth=3, label='Mean trend')
        
        # Plot error bars (standard deviation)
        ax4.fill_between(unique_aspects_sorted, 
                         mean_eff_sorted - std_eff_sorted, 
                         mean_eff_sorted + std_eff_sorted, 
                         alpha=0.2, color='red', label='±1 std dev')

    ax4.set_xlabel('Aspect Ratio (Width/Length)', fontsize=12)
    ax4.set_ylabel('0th Order Diffraction Efficiency (%)', fontsize=12)
    ax4.set_title('Diffraction Efficiency vs Aspect Ratio', fontsize=14, fontweight='bold')
    ax4.tick_params(axis='both', which='major', labelsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=10)
    ax4.set_ylim([0, 100])
else:
    ax4.text(0.5, 0.5, 'No valid data to plot', 
             horizontalalignment='center', verticalalignment='center',
             transform=ax4.transAxes, fontsize=14)
    ax4.set_title('Diffraction Efficiency vs Aspect Ratio', fontsize=14, fontweight='bold')

# Save Plot 4
plot4_path = 'C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\efficiency_vs_aspect_ratio.png'
fig4.savefig(plot4_path, dpi=300, bbox_inches='tight')
print(f">> Saved Plot 4: {plot4_path}")
plt.close(fig4)  # Close the figure to free memory

# Show summary plots (optional)
fig_summary, axes_summary = plt.subplots(2, 2, figsize=(14, 12))

# Recreate small versions for summary
axes_summary[0, 0].imshow(phase_mat, extent=[length_vec[0]*1e9, length_vec[-1]*1e9, 
                                            width_vec[0]*1e9, width_vec[-1]*1e9],
                         aspect='auto', cmap='rainbow', origin='lower')
axes_summary[0, 0].set_title('Phase Shift (rad)')
axes_summary[0, 0].set_xlabel('Length (nm)')
axes_summary[0, 0].set_ylabel('Width (nm)')

axes_summary[0, 1].imshow(diff_eff_mat_plot * 100, extent=[length_vec[0]*1e9, length_vec[-1]*1e9, 
                                                          width_vec[0]*1e9, width_vec[-1]*1e9],
                         aspect='auto', cmap='viridis', origin='lower')
axes_summary[0, 1].set_title('Diffraction Efficiency (%)')
axes_summary[0, 1].set_xlabel('Length (nm)')
axes_summary[0, 1].set_ylabel('Width (nm)')

if len(aspect_flat_valid) > 0:
    axes_summary[1, 0].scatter(aspect_flat_valid, phase_flat_valid, c=diff_eff_flat_valid * 100, 
                              cmap='viridis', s=30, alpha=0.6)
    axes_summary[1, 0].set_title('Phase vs Aspect Ratio')
    axes_summary[1, 0].set_xlabel('Aspect Ratio')
    axes_summary[1, 0].set_ylabel('Phase (rad)')
    axes_summary[1, 0].grid(True, alpha=0.3)

    if len(unique_aspects_sorted) > 0:
        axes_summary[1, 1].plot(unique_aspects_sorted, mean_eff_sorted, 'r-', linewidth=2)
        axes_summary[1, 1].fill_between(unique_aspects_sorted, 
                                       mean_eff_sorted - std_eff_sorted, 
                                       mean_eff_sorted + std_eff_sorted, 
                                       alpha=0.2, color='red')
        axes_summary[1, 1].set_title('Efficiency vs Aspect Ratio')
        axes_summary[1, 1].set_xlabel('Aspect Ratio')
        axes_summary[1, 1].set_ylabel('Efficiency (%)')
        axes_summary[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
summary_path = 'C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\summary_plots.png'
plt.savefig(summary_path, dpi=300, bbox_inches='tight')
plt.show()
plt.close()

# Save all data
print(">> Progress: Saving all data...")
save_dir = 'C:\\Users\\HARSH\\Desktop\\Python_Lumerical\\'
data_path = os.path.join(save_dir, 'rectangular_si_meta_atom_data.npy')

with open(data_path, 'wb') as f:
    np.save(f, width_vec)
    np.save(f, length_vec)
    np.save(f, phase_mat)
    np.save(f, diff_eff_mat)
    np.save(f, aspect_ratio_mat)

# Save as CSV for easier analysis
data_list = []
for i, width in enumerate(width_vec):
    for j, length in enumerate(length_vec):
        if not np.isnan(phase_mat[i, j]):  # Skip failed simulations
            data_list.append({
                'width_nm': width * 1e9,
                'length_nm': length * 1e9,
                'aspect_ratio': width / length,
                'phase_rad': phase_mat[i, j],
                'diffraction_efficiency': diff_eff_mat[i, j],
                'diffraction_efficiency_percent': diff_eff_mat[i, j] * 100
            })

df_results = pd.DataFrame(data_list)
csv_path = os.path.join(save_dir, 'rectangular_si_meta_atom_data.csv')
df_results.to_csv(csv_path, index=False)

print(">> Progress: Simulation completed successfully!")
print("\n=== Summary ===")
print(f"Structure: SiO2 substrate with rectangular Si nano-pillars")
print(f"Si height: {height_atom*1e9:.0f} nm")
print(f"Total simulations: {len(width_vec) * len(length_vec)}")
print(f"Phase range: {np.nanmin(phase_mat):.3f} to {np.nanmax(phase_mat):.3f} rad")
print(f"0th Order Diffraction Efficiency range: {np.nanmin(diff_eff_mat)*100:.1f}% to {np.nanmax(diff_eff_mat)*100:.1f}%")
print(f"\nData saved to:")
print(f"  - {data_path} (NumPy format)")
print(f"  - {csv_path} (CSV format)")
print(f"\nPlots saved to:")
print(f"  1. {plot1_path}")
print(f"  2. {plot2_path}")
print(f"  3. {plot3_path}")
print(f"  4. {plot4_path}")
print(f"  5. {summary_path} (Summary of all plots)")

# Close FDTD
fdtd.close()
print(">> Progress: FDTD closed.")
