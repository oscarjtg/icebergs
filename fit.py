"""
Python functions to fit analytic functions to the ODE solutions 
from the icebergs module.

In development.
"""

import numpy as np
import matplotlib.pyplot as plt
import icebergs

def iceberg_trajectory(height, aspect_ratio, initial_angle_degrees=-5., density_ratio=917./1025., gravity=9.81, n_timesteps=6000):
    """
    Uses iceberg module to calculate the trajectory of a rectangular iceberg 
    which is initially leaning against a smooth vertical wall.

    Scales the time step and drag coefficients appropriately.

    The time step is scaled by sqrt(height / gprime), 
    where 
    gprime = (1 - density_ratio) * gravity 
    is the reduced gravity.

    The x and z drag coefficients are scaled by density_ratio / sqrt(height).

    The angular drag coefficient is scaled by density_ratio.

    Parameters
    ----------
    height (float)
        The height of the iceberg, in m. This is the characteristic lengthscale.

    aspect_ratio (float)
        The ratio of width to height for the iceberg (dimensionless).

    initial_angle_degrees (float)
        The initial angle that the iceberg side of length "height" makes to the vertical, measured clockwise in degrees.
        Default value: -5.0

    density_ratio (float)
        The ratio of densities between the iceberg and the ocean water (dimensionless).
        Default value: 917.0/1025.0

    gravity (float)
        The gravitational acceleration, in m/s^2.
        Default value: 9.81

    n_timesteps (int)
        Integer number of timesteps for which to run the model.

    Returns
    --------
    solver (icebergs.DynamicsSolver)
        Contains 1D numpy arrays of trajectory time series data.
    """
    # Initialise the iceberg shape (`rect`), and the iceberg (`berg`) itself.
    width = aspect_ratio * height
    rect = icebergs.Rectangle(width, height)
    berg = icebergs.Iceberg2D(rect)

    # Set the initial position and orientation of the iceberg.
    theta_deg = initial_angle_degrees
    theta = theta_deg * np.pi / 180
    xwall = 0.0       # The horizontal location of the vertical wall
    x0 = xwall + 0.5 * (height * np.sin(np.abs(theta)) + width * np.cos(theta))

    # Set the density of the iceberg.
    # gprime = g*(rho_w - rho_i)/rho_w
    # therefore rho_i = rho_w * (1 - (gprime / g))
    density_ratio0 = berg.density_water / berg.density_ice
    gprime0 = berg.gravity * (berg.density_water - berg.density_ice) / berg.density_water
    berg.set_density_ice(berg.density_water * density_ratio)
    gprime = (1 - density_ratio) * gravity

    # Set berg.z such that iceberg is in hydrostatic balance.
    berg.set_hydrostatic_balance(x0, theta)

    # Set model parameters.
    # Numerical values for drag coefficients are based on Burton et al (2012).
    height0 = 0.1
    timestep = 0.002 * np.sqrt((height * gprime0) / (height0 * gprime))  # The model time step
    gamma_x = 500. * (height0 / height) * (density_ratio / density_ratio0)    # The horizontal quadratic drag coefficient
    gamma_z = 500. * (height0 / height) * (density_ratio / density_ratio0)    # The vertical quadratic drag coefficient
    gamma_theta = 3. * (density_ratio / density_ratio0)   # The angular quadratic drag coefficient
    restitution = 0.0  # The coefficient of restitution for iceberg-wall collisions
    
    # Configure the solver.
    # n_chunks * n_steps gives the total number of timesteps
    solver = icebergs.DynamicsSolver(timestep, 
                                     gamma_x, 
                                     gamma_z, 
                                     gamma_theta, 
                                     xwall, 
                                     restitution, 
                                     n_timesteps
                                    )

    # Run the solver.   
    solver.simulate(berg, n_timesteps, 
                    plot=False, 
                    saveplot=False
                    )
        
    return solver

def make_dimensionless(t: np.ndarray, 
                       x: np.ndarray, 
                       z: np.ndarray, 
                       theta: np.ndarray, 
                       u: np.ndarray, 
                       w: np.ndarray, 
                       omega: np.ndarray, 
                       height: float, 
                       gprime: float):
    """
    Converts input arrays into arrays of dimensionless quantities.
    
    Parameters
    ----------
    t (np.ndarray)
        1D numpy array contaning times, in seconds.

    x (np.ndarray)
        1D numpy array containing time series data of horizontal coordinate of centre of mass, in metres.

    z (np.ndarray)
        1D numpy array containing time series data of vertical coordinate of centre of mass, in metres.

    theta (np.ndarray)
        1D numpy array contaning time series data of angle, measured clockwise in degrees.

    u (np.ndarray)
        1D numpy array containing time series data of horizontal velocity, in metres.

    w (np.ndarray)
        1D numpy array containing time series data of vertical velocity, in metres.

    omega (np.ndarray)
        1D numpy array containing time series data of angular velocity, in metres.

    height (np.ndarray)
        Height of iceberg.

    gprime (np.ndarray)
        Reduced gravity.

    Returns
    --------
    Tuple of 1D numpy arrays containing non-dimensionlised data from the input arrays.

    """
    t_d = np.sqrt(gprime / height) * t
    x_d = x / height
    z_d = z / height
    theta_d = theta
    u_d = u / np.sqrt(gprime * height)
    w_d = w / np.sqrt(gprime * height)
    omega_d = np.sqrt(height / gprime) * omega
    return t_d, x_d, z_d, theta_d, u_d, w_d, omega_d

def oscillating_decay_model(t, a0, a1, a2, a3):
    """
    Model of the form of the solution of a simple harmonic oscillator 
    with quadratic damping in the weak damping limit (WKB approximation),

    y(t) = a0 / (1 + a1*t) * sin(a2*t + a3)

    Has four parameters for fitting.

    Parameters
    ----------
    t (float)
        The independent variable, time.

    a0 (float)
        Amplitude parameter.

    a1 (float)
        Time amplitude scaling parameter. a1 > 0.

    a2 (float)
        Oscillation frequency parameter.

    a3 (float)
        Oscillation phase parameter.

    Returns
    -------
    y (float)
        The value predicted by the model.
    
    """
    return a0 / (1 + a1 * t) * np.sin(a2 * t + a3)

if __name__ == "__main__":
    height = 1.0
    aspect_ratio = 0.5
    solution = iceberg_trajectory(height, aspect_ratio)
    solution.plot_trajectory()
    plt.show()