import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define the parameters
M = 1
r_0 = 10
u_dot_0 = 0
b = 60
#b = np.sqrt(27)*M

# Define the differential equation
def f(u, phi):
    return np.sqrt(np.abs(2 * M * u**3 - u**2 + 1/b**2))

# Fourth-order Runge-Kutta method
def runge_kutta(h, u, phi):
    k1 = h * f(u, phi)
    k2 = h * f(u + k1/2, phi + h/2)
    k3 = h * f(u + k2/2, phi + h/2)
    k4 = h * f(u + k3, phi + h)
    
    u_new = u + (k1 + 2*k2 + 2*k3 + k4)/6
    phi_new = phi + h
    
    return u_new, phi_new

# Initial conditions
u_0 = 1/r_0
phi_0 = 0

# Time step and number of steps
h = 0.01
num_steps = 10000

# Arrays to store the results
u_values = np.zeros(num_steps)
phi_values = np.zeros(num_steps)

# Perform the integration using Runge-Kutta
u_values[0] = u_0
phi_values[0] = phi_0


for i in range(1, num_steps):
    u_values[i], phi_values[i] = runge_kutta(h, u_values[i-1], phi_values[i-1])

print(u_values)
print(phi_values)

# Convert to Cartesian coordinates
x_values = 1/u_values * np.cos(phi_values)
y_values = 1/u_values * np.sin(phi_values)

fig, ax = plt.subplots()
line, = ax.plot([], [], marker='o')

def init():
    ax.set_xlim(min(x_values), max(x_values))
    ax.set_ylim(min(y_values), max(y_values))
    return line,

def update(frame):
    x_data = x_values[:frame+1]
    y_data = y_values[:frame+1]
    line.set_data(x_data, y_data)
    return line,

# Calcula el número total de frames para 10 segundos a una tasa de 30 fps
total_frames = num_steps
ani = FuncAnimation(fig, update, frames=total_frames, init_func=init, blit=True, interval=10)

plt.show()
