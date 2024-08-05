import numpy as np
import matplotlib.pyplot as plt

# Given parameters
M = 0.8
l = 4.0
r_0 = 35.0
drdt_0 = 0.0

# Define the system of differential equations
def equations(tao, y):
    r, drdt, phi, dphidtao, dt_dtao = y
    
    # Define the derivatives
    drdtao = drdt
    ddrdt_dtao = -M / r**2 + l**2 / r**3 - 2 * M * l**2 / r**4
    dphidtao = l / r**2
    ddt_dtao = np.sqrt(r / (r - 2 * M) + r**2 * dphidtao**2 / (r - 2 * M) + r**2 * drdt**2 / (r - 2)**2)
    
    return [drdt, ddrdt_dtao, dphidtao, 0, ddt_dtao]

# Implement the fourth-order Runge-Kutta method
def runge_kutta(h, tao_max):
    num_steps = int(tao_max / h)
    tao_values = np.linspace(0, tao_max, num_steps + 1)
    y_values = np.zeros((num_steps + 1, 5))  # 5 variables: r, drdt, phi, dphidtao, dt_dtao
    
    # Initial conditions
    y_values[0] = [r_0, drdt_0, 0, 0, np.sqrt(r_0 / (r_0 - 2 * M))]
    
    for i in range(num_steps):
        k1 = h * np.array(equations(tao_values[i], y_values[i]))
        k2 = h * np.array(equations(tao_values[i] + h / 2, y_values[i] + k1 / 2))
        k3 = h * np.array(equations(tao_values[i] + h / 2, y_values[i] + k2 / 2))
        k4 = h * np.array(equations(tao_values[i] + h, y_values[i] + k3))
        
        y_values[i + 1] = y_values[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
    
    return tao_values, y_values

# Set the time step and maximum tao value
h = 0.1
tao_max = 9000

# Run the simulation
tao_values, y_values = runge_kutta(h, tao_max)

# Convert polar coordinates to Cartesian coordinates
x_values = y_values[:, 0] * np.cos(y_values[:, 2])
y_values = y_values[:, 0] * np.sin(y_values[:, 2])

# Plot the results in xy plane
plt.figure(figsize=(8, 8))
plt.plot(x_values, y_values, label='Orbit in xy plane')
plt.title('Orbit in xy plane')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()
