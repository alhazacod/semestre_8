import numpy as np
import matplotlib.pyplot as plt

# Define the system of equations
def equations(r, phi, m, b):
    dr_dphi = np.sqrt(-2 * m * r + r**2 - r**4 / b**2)
    return dr_dphi

# Define the fourth-order Runge-Kutta method
def runge_kutta(h, r, phi, m, b):
    k1 = h * equations(r, phi, m, b)
    k2 = h * equations(r + 0.5 * k1, phi + 0.5 * h, m, b)
    k3 = h * equations(r + 0.5 * k2, phi + 0.5 * h, m, b)
    k4 = h * equations(r + k3, phi + h, m, b)
    return r + (k1 + 2 * k2 + 2 * k3 + k4) / 6

# Initial conditions
r_0 = 5
phi_0 = 0
m = 1
b = 4

# Time step and number of steps
h = 0.01
num_steps = 1000

# Arrays to store the results
r_values = np.zeros(num_steps)
phi_values = np.zeros(num_steps)

# Perform the Runge-Kutta integration
for i in range(num_steps):
    r_values[i] = r_0
    phi_values[i] = phi_0
    r_0 = runge_kutta(h, r_0, phi_0, m, b)
    phi_0 += h

# Convert to Cartesian coordinates
x_values = r_values * np.cos(phi_values)
y_values = r_values * np.sin(phi_values)

# Plot the solution
plt.plot(x_values, y_values)
plt.title('Solution of the System of Equations')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

