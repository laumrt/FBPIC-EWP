import math
import numpy as np
from numba import jit, prange
from scipy.interpolate import interp1d
from os import environ
environ['NUMBA_NUM_THREADS'] = '8'

@jit(nopython=True)
def interp_2D(x, y, xmin, ymin, dx, dy, E1, E2):
    """
    Perform 2D bilinear interpolation to interpolate a value from a 2D grid.

    Parameters:
    - x, y: Coordinates of the point to be interpolated.
    - xmin, ymin: Minimum values of x and y in the grid.
    - dx, dy: Grid spacings in the x and y directions.
    - E1, E2: 2D arrays representing the grid of values to interpolate from.

    Returns:
    - E1_interp, E2_interp: Interpolated values at the given coordinates (x, y).
    """
    Nx, Ny = E1.shape
    assert E1.shape==E2.shape

    # Calculate indices and coefficients in x direction
    s_ix = (x - xmin) / dx   # Calculate fractional index in x direction
    ix = np.floor(s_ix)      # Integer part of the index
    sx1 = s_ix - ix          # Fractional part of the index
    sx0 = 1.0 - sx1          # Complementary fraction
    ix = int(ix)             # Convert to integer

    # Calculate indices and coefficients in y direction
    s_iy = (y - ymin) / dy   # Calculate fractional index in y direction
    iy = np.floor(s_iy)      # Integer part of the index
    sy1 = s_iy - iy          # Fractional part of the index
    sy0 = 1.0 - sy1          # Complementary fraction
    iy = int(iy)             # Convert to integer

    # Perform boundary checks
    if ix<0 or ix>=Nx-1 or iy<0 or iy>=Ny-1:
        return 0.0, 0.0  # Return 0 if indices are out of bounds

    # Calculate the weights for bilinear interpolation
    C00 = sx0 * sy0   # Weight for the lower left corner
    C01 = sx0 * sy1   # Weight for the upper left corner
    C10 = sx1 * sy0   # Weight for the lower right corner
    C11 = sx1 * sy1   # Weight for the upper right corner

    # Perform bilinear interpolation
    E1_interp = (E1[ix, iy] * C00 +
         E1[ix, iy + 1] * C01 +
         E1[ix + 1, iy] * C10 +
         E1[ix + 1, iy + 1] * C11
    )

    # Perform bilinear interpolation
    E2_interp = (E2[ix, iy] * C00 +
         E2[ix, iy + 1] * C01 +
         E2[ix + 1, iy] * C10 +
         E2[ix + 1, iy + 1] * C11
    )

    return E1_interp, E2_interp

@jit(nopython=True, parallel=True)
def particles_loop(x, y, z, Ex_t, Ey_t, Ez_t, r_int_min, z_int_min, dr_int, dz_int, Er_w_map, Ez_w_map) :
    """
    Perform particle loop calculations to interpolate and compute electric field components.

    Parameters:
    - x, y, z: Arrays of particle coordinates in Cartesian coordinates (m).
    - Ex_t, Ey_t, Ez_t: Arrays of total electric field components at particle positions (V/m).

    Returns:
    - Ex_l, Ey_l, Ez_l: Arrays of laser-induced electric field components (V/m).
    - Ex_w, Ey_w, Ez_w: Arrays of wake electric field components (V/m). Each array entry corresponds to the electric field on a particle at a certain iteration.
    """
    N_p = x.size
    Ex_w, Ey_w, Ez_w = np.zeros(N_p), np.zeros(N_p), np.zeros(N_p)
    Ex_l, Ey_l, Ez_l = np.zeros(N_p), np.zeros(N_p), np.zeros(N_p)

    for ip in prange(N_p) :
        x_p, y_p, z_p, Ex_t_p, Ey_t_p, Ez_t_p = x[ip], y[ip], z[ip], Ex_t[ip], Ey_t[ip], Ez_t[ip]

        # Calculate radial distance and angle
        r_p = math.sqrt(x_p**2 + y_p**2) # m
        theta_p = math.atan2(y_p, x_p) # rad

        # Interpolate electric field components at particle position (V/m)
        Er_w_int, Ez_w_int = interp_2D(
            r_p, z_p, r_int_min, z_int_min,
            dr_int, dz_int, Er_w_map, Ez_w_map
        )

        # Transform cylindrical electric field components to Cartesian coordinates (V/m)
        Ex_w[ip] = Er_w_int * np.cos(theta_p)
        Ey_w[ip] = Er_w_int * np.sin(theta_p)
        Ez_w[ip] = Ez_w_int

        # Calculate laser-induced electric field components (V/m)
        Ex_l[ip] = Ex_t_p - Ex_w[ip]
        Ey_l[ip] = Ey_t_p - Ey_w[ip]
        Ez_l[ip] = Ez_t_p - Ez_w[ip]

    return Ex_l, Ey_l, Ez_l, Ex_w, Ey_w, Ez_w

def spatial_lowpass_filter(field, dz, lambda0=0.8e-6, cutoff_freq=0.5):
    """
    Apply a low-pass filter to the 1D spatial field along the x-direction using the Fourier transform.

    Parameters:
    - field: 1D array representing the spatial field to be filtered.
    - cutoff_freq: Normalized Cutoff frequency for the low-pass filter, specified in terms of spatial frequencies. Default is 0.5.
    - dz: Spatial sampling interval of the field along the x-direction.
    - lambda0: Wavelength corresponding to the central frequency of the field, in meters. Default is 0.8e-6 meters.

    Returns:
    - filtered_field: 1D array representing the filtered spatial field. The field is filtered to attenuate spatial frequencies higher than the cutoff frequency while preserving lower spatial frequencies.
    """
    # Compute the maximum value of the input field
    max_value = np.max(field)

    k0 = 2 * np.pi / lambda0

    # Compute the 1D Fourier transform of the spatial field
    field_fft = np.fft.rfft(field)

    # Compute the frequency axis
    freq = 2 * np.pi * np.fft.rfftfreq(field.shape[1], d=dz) / k0

    # Apply the low-pass filter in the frequency domain
    field_fft *= (freq<cutoff_freq)[None, :]

    # Compute the inverse Fourier transform to obtain the filtered spatial field
    filtered_field = np.real(np.fft.irfft(field_fft))

    # Scale the filtered field to have the same maximum value as the input field
    filtered_field *= max_value / np.max(filtered_field)

    return filtered_field

def shift_half_step(A, t, dt):
    """
    Shifts the input array `A` by half a time step using cubic interpolation.

    Parameters:
    - A: 1D array_like. The input array to be shifted.
    - t: 1D array_like. The time values corresponding to the data in `A`.
    - dt: float. The time step interval.

    Returns:
    - A_new: 1D array. The shifted array after performing cubic interpolation.

    Notes:
    This function performs cubic interpolation on the input array `A` using the time values `t`.
    It then shifts the interpolated values forward by half a time step (`0.5 * dt`) and returns
    the interpolated array at the shifted time points.
    """
    interp_fu = interp1d(t, A, assume_sorted=True,
                         kind='cubic', bounds_error=None,
                         fill_value='extrapolate')
    A_new = interp_fu(t + 0.5 * dt)
    return A_new

