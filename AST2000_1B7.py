# IKKE KODEMAL
# SKREVET EGEN KODE :)
# (Python 3.9)

# Import neccesary packages
from time import time 
start_time = time() # The time of when the code starts running

# Packages used for computations
import numpy as np

# Packages from ast2000tools that we need
import ast2000tools.utils as ast
import ast2000tools.constants as const
from ast2000tools.solar_system import SolarSystem

# Packages for visualization of the simulation
import matplotlib.pyplot as plt

# Constants
G:float = const.G # Gravtational constant
m_sun:float = const.m_sun # Solar mass
AU:float = const.AU # Astronomical unit
year:int = 365.2421988*24*60*60 #Seconds in a year

# Creating a random solar system with AST2000tools
seed:int = ast.get_seed('abdulaah')
solarsystem:SolarSystem = SolarSystem(seed)

#Initial values
frames:int = 1000 # Amount of steps 
deltaT:float = 1000000 # Length of each step

star_m:float = solarsystem.star_mass*m_sun # Star mass
planet_m:float = solarsystem.masses*m_sun # Orbiting planet masses
initial_positions:float = solarsystem.initial_positions*AU # Initial position of planets in meters
initial_velocitites:float = solarsystem.initial_velocities*(AU/year) # Initial velocites of planet in m/s

print(initial_positions)
print(initial_velocitites)

planet_positions:np.ndarray = np.zeros((frames, np.shape(initial_positions)[0], np.shape(initial_positions)[1])) # Array for x and y positions of each planet
planet_velocities:np.ndarray = np.zeros((frames, np.shape(initial_velocitites)[0], np.shape(initial_velocitites)[1])) # Array for Vx and Vy velocities of each planet
planet_positions[0] = initial_positions # Setting initial positions 
planet_velocities[0] = initial_velocitites # Setting initial velocities

# Using Euler-Cromer method to find planets orbit for each timestep
for i in range(frames - 1):
    for j in range(len(planet_m)):
        x = planet_positions[i, 0, j]
        y = planet_positions[i, 1, j]
        r = np.sqrt(x**2 + y**2) # Distance between star and planet
        g = (G * star_m) / (r**2) # Gravitational akseleration between planet and star

        # Calculationg aksleration in x and y
        gx = -g*(x/r)
        gy = -g*(y/r)

        # Updating velocities and positions for the next frame
        planet_velocities[i + 1, 0, j] = planet_velocities[i, 0, j] + gx*deltaT
        planet_velocities[i + 1, 1, j] = planet_velocities[i, 1, j] + gy*deltaT

        planet_positions[i + 1, 0, j] = planet_positions[i, 0, j] + planet_velocities[i + 1, 0, j]*deltaT
        planet_positions[i + 1, 1, j] = planet_positions[i, 1, j] + planet_velocities[i + 1, 1, j]*deltaT

# Calculating potential and kinetic energies of the planets
potential_energies:np.ndarray = -(G*star_m*planet_m)/(np.sqrt(planet_positions[:, 0, :]**2 + planet_positions[:, 1, :]**2))
kinetic_energies:np.ndarray = (1/2)*planet_m*(planet_velocities[:, 0, :]**2 + planet_velocities[:, 1, :]**2)
total_energies = kinetic_energies + potential_energies

end_time = time() # The time of when the code stops running
print()
print(f'The code used \033[1m{end_time - start_time}\033[0m seconds to fully execute')
print()

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink']
planet = 2

# Plotting the orbits
'''
#plt.style.use('dark_background')
plt.figure(figsize=(8, 8))
for j in range(len(planet_m)):
    plt.plot(planet_positions[:, 0, j], planet_positions[:, 1, j], color=colors[j], label=f'Planet {j + 1} (Mass: {planet_m[j]:.2e} kg)')
    plt.arrow(planet_positions[frames - 2, 0, j], planet_positions[frames - 2, 1, j], (planet_positions[frames - 1, 0, j] - planet_positions[frames - 2, 0, j]), (planet_positions[frames - 1, 1, j] - planet_positions[frames - 2, 1, j]), head_width=1e11, head_length=1e11, fc=colors[j], ec=colors[j])

# Plot the star at the center
plt.scatter(0, 0, color='yellow', s=200, label=f'Star       (Mass: {star_m:.2e} kg)')

# Labels and legend
plt.xlabel('x-position [m]')
plt.ylabel('y-position [m]')
plt.title('Planetary Orbits')
plt.legend()
plt.grid()
plt.show()
'''


# Plotting energies
'''
plt.figure(figsize=(7, 10))

# Potential Energy
plt.subplot(3, 1, 1)
plt.title(f'Planet {planet} (delta T = 10 000)')
plt.plot(np.linspace(0, frames*deltaT, frames), potential_energies[:,planet - 1], color=colors[planet - 1])
plt.ylabel('Potential Energy (J)')
plt.grid()

# Kinetic Energy
plt.subplot(3, 1, 2)
plt.plot(np.linspace(0, frames*deltaT, frames), kinetic_energies[:,planet - 1], color=colors[planet - 1])
plt.ylabel('Kinetic Energy (J)')
plt.grid()

# Total Energy
centre = (min(total_energies[:,planet - 1]) + max(total_energies[:,planet - 1]))/2
scaling = (max(potential_energies[:,planet - 1]) - min(potential_energies[:,planet - 1]) + max(kinetic_energies[:,planet - 1]) - min(kinetic_energies[:,planet - 1]))/4

plt.subplot(3, 1, 3)
plt.plot(np.linspace(0, frames*deltaT, frames), total_energies[:,planet - 1], color=colors[planet - 1])
plt.xlabel('Time (s)')
plt.ylabel('Total Energy (J)')
plt.grid()
plt.ylim(centre - scaling, centre + scaling)

plt.show()
'''

# Plotting change in radius
'''
plt.figure(figsize=(7, 7))
plt.title(f'Planet {planet}')
plt.plot(np.linspace(0, frames*deltaT, frames), np.sqrt(planet_positions[:, 0, planet - 1]**2 + planet_positions[:, 1, planet - 1]**2), color=colors[planet - 1])
plt.ylabel('Distance between Planet and Star (m)')
plt.xlabel('Time (s)')
plt.grid()

plt.show()
'''