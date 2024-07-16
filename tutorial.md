# Code Explanation (main.py)
This tutorial guides you through the process of analyzing FBPIC simulation results simulation data using FBPIC-EWP. Here, we provide an explanation of the different sections of the code, allowing the user to define the parameters for the postprocessing analysis.

## 1. Imports and Initial Setup
First, we import the necessary libraries and define the output file name (`file_out`) and simulation path (`hdf5_path`). Subsequently, the code imports the simulation (`sim = LpaDiagnostics(...)`) and it is possible to set `magnetic_field = True`, in case to calculate the laser and plasma magnetic field on each tracked particle. 

```python
import numpy as np
from scipy.constants import c, e, m_e, pi
import time
import h5py
from tqdm import tqdm
from openpmd_viewer import OpenPMDTimeSeries, ParticleTracker
from openpmd_viewer.addons import LpaDiagnostics
from func import particles_loop, spatial_lowpass_filter, shift_half_step

# Defining the output file name
file_out = 'work_result_test.h5'

# Selecting the simulation
hdf5_path = ".../hdf5"
sim = LpaDiagnostics(hdf5_path, backend='h5py', check_all_files=False)

# Set to True if you want to decompose the magnetic field in m=0,1 modes (False by default)
magnetic_field = False
```
## 2. Particle Selection
Select the species of particles to study that has to correspond to a species defined in the original FBPIC simulation (in this example `species='elec_L'`). Here, we import the ids of the tracked particles at the last iteration with `id_ = sim.get_particle(['id'], t = sim.t[-1], species=species)[0]`. We then randomly choose a subset of ten particles to study (i.e., `id_array = np.random.choice(id_, size=10, replace=False)`).
```python
# Selecting the particles to study
species = 'elec_L'
id_ = sim.get_particle(['id'], t = sim.t[-1], species=species)[0]
id_array = np.random.choice(id_, size=10, replace=False)
pt_select = ParticleTracker(sim, species=species, preserve_particle_index=True, select=id_array)
N_p = pt_select.N_selected  # number of selected electrons
```
## 3. Timestep Parameters
FBPIC-EWP automatically defines the timestep (`dt`) from the original FBPIC simulation timestep.
```python
# Define the time range
t = sim.t  # (s)
N_t = t.size
dt = t[0:2].ptp()  # (s), time step
```
## 4. Data Import
Import the coordinates (`x, y, z`), momenta (`ux, uy, uz`), and total electric (`Ex_t, Ey_t, Ez_t`) and magnetic (`Bx_t, By_t, Bz_t`) fields of each particle for every iteration.
```python
# Importing the coordinates, momenta and total EM-field of each particle for every iteration
if magnetic_field:
    print('(1/3) Importing the coordinates, momenta and total EM-field of each particle for every iteration\n')
    x, y, z, ux, uy, uz, Ex_t, Ey_t, Ez_t, Bx_t, By_t, Bz_t, w = sim.iterate(sim.get_particle,
       ['x', 'y', 'z', 'ux', 'uy', 'uz', 'E/x', 'E/y', 'E/z', 'B/x', 'B/y', 'B/z', 'w'],
       spec=species, select=pt_select)
    Bx_t, By_t, Bz_t = np.nan_to_num(Bx_t), np.nan_to_num(By_t), np.nan_to_num(Bz_t)
else:
    print('(1/3) Importing the coordinates, momenta and total E-field of each particle for every iteration\n')
    x, y, z, ux, uy, uz, Ex_t, Ey_t, Ez_t, w = sim.iterate(sim.get_particle,
       ['x', 'y', 'z', 'ux', 'uy', 'uz', 'E/x', 'E/y', 'E/z', 'w'],
       spec=species, select=pt_select)
```
## 5. Field Calculation Loop
For each time step, extract and filter the electric and magnetic field components. Use a JIT-compiled loop (see `func.py`) to interpolate the fields onto the particle positions and save the electric field components for each particle. Here, (`Ex_las, Ey_las, Ez_las`) and (`Ex_wake, Ey_wake, Ez_wake`) are arrays containing the laser and plasma electric field at each iteration for each selected particle. Analogously, (`Bx_las, By_las, Bz_las`) and (`Bx_wake, By_wake, Bz_wake`) are arrays containing the laser and plasma magnetic field at each iteration for each selected particle. 
```python
print('(2/3) Calculating the electric fields on each particle... \n')

# Electric field on particle calculation (iterate over each time step in the specified time range)
for it in tqdm(range(N_t)):
    
    ################################ ELECTRIC FIELD ##########################################
    
    # Extract electric field components in cylindrical coordinates (V/m)
    Er_w_map = sim.get_field(t=t[it], field='E', coord='r', m=0, theta=pi/2)[0]
    Ez_w_map, info = sim.get_field(t=t[it], field='E', coord='z', m=0, theta=pi/2)

    # Lowpass filter on the radial field
    Er_w_map = spatial_lowpass_filter(Er_w_map, info.dz)

    # Extract radial and axial coordinates for field interpolation
    r_int, dr_int = info.r, info.dr  # m
    z_int, dz_int = info.z, info.dz  # m
    
    # JIT particle loop
    Ex_l, Ey_l, Ez_l, Ex_w, Ey_w, Ez_w = particles_loop(
        x[it, :], y[it, :], z[it, :], Ex_t[it, :], Ey_t[it, :], Ez_t[it, :], 
        r_int.min(), z_int.min(), dr_int, dz_int,
        Er_w_map, Ez_w_map
    )

    # Saving all the electric fields (V/m)
    Ex_las[it, :] = Ex_l
    Ey_las[it, :] = Ey_l
    Ez_las[it, :] = Ez_l
    Ex_wake[it, :] = Ex_w
    Ey_wake[it, :] = Ey_w
    Ez_wake[it, :] = Ez_w
    
    ################################# MAGNETIC FIELD ##########################################
    
    if magnetic_field:
        # Extract magnetic field components in cylindrical coordinates (T)
        Br_w_map = sim.get_field(t=t[it], field='B', coord='r', m=0, theta=pi/2)[0]
        Bz_w_map, info = sim.get_field(t=t[it], field='B', coord='z', m=0, theta=pi/2)

        # Lowpass filter on the radial field
        Br_w_map = spatial_lowpass_filter(Br_w_map, info.dz)

        # Extract radial and axial coordinates for field interpolation
        r_int, dr_int = info.r, info.dr  # m
        z_int, dz_int = info.z, info.dz  # m

        # JIT particle loop
        Bx_l, By_l, Bz_l, Bx_w, By_w, Bz_w = particles_loop(
            x[it, :], y[it, :], z[it, :], Bx_t[it, :], By_t[it, :], Bz_t[it, :],
            r_int.min(), z_int.min(), dr_int, dz_int,
            Br_w_map, Bz_w_map
        )

        # Saving all the magnetic fields (T)
        Bx_las[it, :] = Bx_l
        By_las[it, :] = By_l
        Bz_las[it, :] = Bz_l
        Bx_wake[it, :] = Bx_w
        By_wake[it, :] = By_w
        Bz_wake[it, :] = Bz_w

```
## 6. Work Calculation
For each particle, shift the momenta by half a time step, calculate the Lorentz factor and normalized velocities, compute the work done by the laser and wakefield components, and save the results. Here, (`Wx_las, Wy_las, Wz_las`) and (`Wx_wake, Wy_wake, Wz_wake`) are arrays containing the laser and plasma work at each iteration for each selected particle.
```python
print('(3/3) Calculating the laser and wakefield work on each particle... \n')

########################## Solution ##########################
for ip in tqdm(range(N_p)):
    ux_p, uy_p, uz_p = ux[:, ip], uy[:, ip], uz[:, ip]
    Ex_las_p, Ey_las_p, Ez_las_p = Ex_las[:, ip], Ey_las[:, ip], Ez_las[:, ip]
    Ex_wake_p, Ey_wake_p, Ez_wake_p = Ex_wake[:, ip], Ey_wake[:, ip], Ez_wake[:, ip]

    ux_p = shift_half_step(ux_p, t, dt)
    uy_p = shift_half_step(uy_p, t, dt)
    uz_p = shift_half_step(uz_p, t, dt)
    
    gamma_p = (1. + ux_p**2 + uy_p**2 + uz_p**2)**0.5
    bx_p, by_p, bz_p = ux_p / gamma_p, uy_p / gamma_p, uz_p / gamma_p

    Wx_las_p = -e * c * np.cumsum(Ex_las_p * bx_p) * dt
    Wy_las_p = -e * c * np.cumsum(Ey_las_p * by_p) * dt
    Wz_las_p = -e * c * np.cumsum(Ez_las_p * bz_p) * dt

    Wx_wake_p = -e * c * np.cumsum(Ex_wake_p * bx_p) * dt
    Wy_wake_p = -e * c * np.cumsum(Ey_wake_p * by_p) * dt
    Wz_wake_p = -e * c * np.cumsum(Ez_wake_p * bz_p) * dt
    
    Wx_las[:, ip] = Wx_las_p
    Wy_las[:, ip] = Wy_las_p
    Wz_las[:, ip] = Wz_las_p

    Wx_wake[:, ip] = Wx_wake_p
    Wy_wake[:, ip] = Wy_wake_p
    Wz_wake[:, ip] = Wz_wake_p
    
    gamma_c[:, ip] = gamma_p
```

## 8. Save Results to HDF5
Finally, save all the calculated data to an HDF5 file.
```python
# Saving the data in a .h5 file
f_out = h5py.File(file_out, 'w')

f_out['w'], f_out['ids'] = w, id_array
f_out['x'], f_out['y'], f_out['z'] = x, y, z
f_out['ux'], f_out['uy'], f_out['uz'] = ux, uy, uz

f_out['Ex_wake'], f_out['Ey_wake'], f_out['Ez_wake'] = Ex_wake, Ey_wake, Ez_wake
f_out['Ex_las'], f_out['Ey_las'], f_out['Ez_las'] = Ex_las, Ey_las, Ez_las

f_out['Wx_wake'], f_out['Wy_wake'], f_out['Wz_wake'] = Wx_wake, Wy_wake, Wz_wake
f_out['Wx_las'], f_out['Wy_las'], f_out['Wz_las'] = Wx_las, Wy_las, Wz_las

if magnetic_field:
    f_out['Bx_wake'], f_out['By_wake'], f_out['Bz_wake'] = Bx_wake, By_wake, Bz_wake
    f_out['Bx_las'], f_out['By_las'], f_out['Bz_las'] = Bx_las, By_las, Bz_las

f_out.close()
```
