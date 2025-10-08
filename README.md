# icebergs
Simple model of floating iceberg dynamics.

## Installation instructions

Clone the repository, place `icebergs.py` in the same directory as your scripts 
(or add its directory to you Python path)
and import the module using
```python
import icebergs
```

## Simple usage example
```python
import icebergs
import numpy as np

def init_depth(rho_i, rho_w, H):
    """Returns an appropriate initial depth for the iceberg, based on hydrostatic balance"""
    return (0.5 - (rho_i/rho_w)) * H

# Initialise the iceberg shape (`rect`), and the iceberg (`berg`) itself.
width = 0.05; height = 0.1
rect = icebergs.Rectangle(width, height)
berg = icebergs.Iceberg2D(rect)

# Set the initial position and orientation of the iceberg.
theta_deg = -1.0
theta = theta_deg * np.pi / 180
z_init = init_depth(berg.density_ice, berg.density_water, rect.b)
berg.set_position(x=0.026, z=z_init, theta=theta)

# Set model parameters
timestep = 0.002   # The model time step
gamma_u = 500.     # The horizontal quadratic drag coefficient
gamma_w = 500.     # The vertical quadratic drag coefficient
gamma_omega = 3.   # The angular quadratic drag coefficient
xwall = 0.0        # The horizontal location of the vertical wall
restitution = 0.0  # The coefficient of restitution for iceberg-wall collisions

# Configure the solver.
# n_chunks * n_steps gives the total number of timesteps
n_chunks = 300
n_steps = 10
solver = icebergs.DynamicsSolver(timestep, 
                                 gamma_u, 
                                 gamma_w, 
                                 gamma_omega, 
                                 xwall, 
                                 restitution, 
                                 n_chunks*n_steps
                                )

# Run the solver.
# If plot=True, displays a plot after every `n_steps` timesteps, producing `n_chunks` plots.
# If saveplot=True, saves these plots as PNG files to the `./plots/iceberg` directory.
# The suggested plotlims makes a nice video for these parameter values.
plotlims = [-width, 1.5 * (width + height), -1.5 * height, 1.5 * height]
for _ in range(n_chunks):    
    solver.simulate(berg, n_steps, 
                    plot=True, 
                    saveplot=True, 
                    plotlims=plotlims
                   )

# Optional: plot the trajectory.
# Produces separate plots of x, z, theta, u, w, and omega against time.
solver.plot_trajectory()
plt.savefig("./plots/trajectory.png", dpi=150)

# Optional, convert the plots to a video.
# The video will be saved at `./videos/iceberg.mp4`
# Required: otgraph, my plotting utilities package, available at https://github.com/oscarjtg/otgraph
import otgraph.video as vid
vid.pngs_to_mp4("./plots/iceberg", "./videos/iceberg", fps=30)
```
