import numpy as np
from matplotlib import pyplot as plt 

# Set the channel dimension constants (in m)
W = 0.001
H = 0.0005
L = 0.008 

# Set the fluid viscosity (in Pa.s)
mu = 1e-3

# Set the height and width rearrangements
w = W/2
h = H/2

# Setting up the transmissability functon
def transmissability(k):
    
    constant = ((H**3.)*W)/(12.*mu)
    shear = []

    for n in range(1,k):
        beta_n = (2*n-1)*((np.pi)/2.0)
        gamma = w/h
        sum_1 = (np.tanh(beta_n*gamma)/beta_n**5.)
        shear.append(sum_1)
        sum_2 = 1. - (6./gamma)*np.sum(shear)

    return constant*sum_2

# Setting up the velocity flow function
def velocity_flow(x, y, Q, k):
    
    constant = (Q/transmissability(k))*((h**2)/(2.*mu))
    shear = []

    for n in range(1,k):
        beta_n = (2*n-1)*((np.pi/2.0))
        gamma = w/h
        sum_1 = (((-1.)**n)/beta_n**2.)*(((np.cosh(beta_n*(x/h)))/(np.cosh(beta_n*gamma)))*np.cos(beta_n*(y/h)))
        shear.append(sum_1)
        sum_2 = (1. - (y/h)**2 + 4*(np.sum(shear)))
    return constant * sum_2

# Setting up the function for max shear stress acting on the top wall
def max_shear_stress_top_wall(k,Q):
    
    constant = (Q/transmissability(k))*H#*((h**2)/(2.*mu))
    shear = []

    for n in range(1,k):
        beta_n = (2*n-1)*((np.pi)/2.0)
        gamma = w/h
        sum_1 = ((-1.)**n/beta_n**2.)*(np.sin(beta_n)/np.cosh(beta_n*gamma))
        shear.append(sum_1)
        sum_2 = (-1./2) - np.sum(shear)
    return constant * sum_2

# Update the x-axis range to 20 ml/day
Q_ml_day = np.linspace(0, 25, 100)
Q = Q_ml_day / (86400. * 1.0e6)  # Convert from ml/day to m^3/s

# Calculate max shear stress for each flow ratecle
shear_value = np.zeros(100)
for i in range(100):
    shear_value[i] = abs(max_shear_stress_top_wall(100, Q[i]))

# Plot the graph
plt.plot(Q_ml_day, shear_value)
plt.xlabel('Flow in ml per day')
plt.ylabel('Shear stress (Pa)')

# Shear stress values for reference
# 1.0mm (W) x 0.5mm (H) x 8.0mm (L)
#At 20.52 ml/day (0.00095)
max_shear_1 = 0.019372113
#At 18.36 ml/day (0.00085)
max_shear_2 = 0.017332587
#At 21.384 ml/day (0.00099)
max_shear_3 = 0.020187948

# Add horizontal lines for mean, max, and min shear stress values
plt.axhline(y=max_shear_1, color='red', linestyle='--', label=f'Max Shear 20.52 ml/day - 1.0mmx0.5mm: {max_shear_1} Pa')
plt.axhline(y=max_shear_2, color='orange', linestyle='-.', label=f'Max Shear 18.36 ml/day - 1.0mmx0.5mm: {max_shear_2} Pa')
plt.axhline(y=max_shear_3, color='green', linestyle='-.', label=f'Max Shear 21.384 ml/day - 1.0mmx0.5mm: {max_shear_3} Pa')

# Add markers at # ml/day
plt.scatter([20.52], [max_shear_1], color='red')
plt.scatter([18.36], [max_shear_2], color='orange')
plt.scatter([21.384], [max_shear_3], color='green')

# Create plot legend
plt.legend()  
# Setting x-axis limit to 25 ml/day
plt.xlim(0, 25)  
# Setting y-axis limit to 0.05 Pa
plt.ylim(0.00, 0.05)
plt.show() 

'''
import numpy as np
from matplotlib import pyplot as plt 

# Set the channel dimension constants (in m)
W = 0.003
H = 0.0005

# Set the fluid viscosity (in Pa.s)
mu = 1e-3

# Set the height and width rearrangements
w = W / 2
h = H / 2

# Function to handle large cosh(x) values by using an approximation
def safe_cosh(x):
    # For large values of x, use approximation with log scaling to avoid overflow
    if np.abs(x) > 700:  # Threshold to prevent overflow in exp
        return np.inf  # cosh(x) approaches infinity for very large x
    else:
        return np.cosh(x)  # Normal computation for smaller values

# Setting up the transmissability function
def transmissability(k):
    constant = ((H**3) * W) / (12 * mu)
    
    for i in range(1, k):
        beta_k = (2 * i - 1) * (np.pi / 2.0)
        gamma = w / h
        sum_1 = (np.tanh(beta_k * gamma) / beta_k**5)
        sum_2 = 1. - (6. / gamma) * np.sum(sum_1)
    
    return constant * sum_2

# Setting up the velocity flow function
def velocity_flow(x, y, Q, k):
    constant = (Q / transmissability(k)) * ((h**2) / (2.*mu))
    
    for i in range(1, k):
        beta_k = (2 * i - 1) * (np.pi / 2.0)
        gamma = w / h
        sum_1 = (((-1.)**i) / beta_k**2) * (safe_cosh(beta_k * (x / h)) / safe_cosh(beta_k * gamma)) * np.cos(beta_k * (y / h))
        sum_2 = (1. - (y / h)**2 + 4 * np.sum(sum_1))
    
    return constant * sum_2

# Setting up the function for max shear stress acting on the top wall
def max_shear_stress_top_wall(k, Q):
    constant = (Q / transmissability(k)) * ((h**2) / (2 * mu))
    
    for i in range(1, k):
        beta_k = (2 * i - 1) * (np.pi / 2.0)
        gamma = w / h
        sum_1 = ((-1.)**i / beta_k**2) * (np.sin(beta_k) / safe_cosh(beta_k * gamma))
        sum_2 = (-1. / 2) - np.sum(sum_1)
    
    return mu * constant * sum_2

# Update the x-axis range to 20 ml/day
Q_ml_day = np.linspace(0, 20, 100)
Q = Q_ml_day / (86400. * 1.0e6)  # Convert from ml/day to m^3/s

# Calculate max shear stress for each flow rate
shear_value = np.zeros(100)
for i in range(100):
    shear_value[i] = abs(max_shear_stress_top_wall(100, Q[i]))

# Plot the graph
plt.plot(Q_ml_day, shear_value)
plt.xlabel('Flow in ml per day')
plt.ylabel('Shear stress (Pa)')

# Shear stress values for reference
# 1.0mm (W) x 0.5mm (H) x 8.0mm (L)
# At 20.52 ml/day (0.00095)
max_shear_1 = 0.019372113
# At 18.36 ml/day (0.00085)
max_shear_2 = 0.017332587
# At 21.384 ml/day (0.00099)
max_shear_3 = 0.020187948

# Add horizontal lines for mean, max, and min shear stress values
plt.axhline(y=max_shear_1, color='red', linestyle='--', label=f'Max Shear 20.52 ml/day - 1.0mmx0.5mm: {max_shear_1} Pa')
plt.axhline(y=max_shear_2, color='orange', linestyle='-.', label=f'Max Shear 18.36 ml/day - 1.0mmx0.5mm: {max_shear_2} Pa')
plt.axhline(y=max_shear_3, color='green', linestyle='-.', label=f'Max Shear 21.384 ml/day - 1.0mmx0.5mm: {max_shear_3} Pa')

# Add markers at the specific flow rates
plt.scatter([20.52], [max_shear_1], color='red')
plt.scatter([18.36], [max_shear_2], color='orange')
plt.scatter([21.384], [max_shear_3], color='green')

# Create plot legend
plt.legend()  
# Setting x-axis limit to 25 ml/day
plt.xlim(0, 25)  
# Setting y-axis limit to 0.05 Pa
plt.ylim(0, 0.05)
plt.show()

'''
