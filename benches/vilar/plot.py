import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

with open('vilar.csv', 'r') as file:
    header = file.readline().rstrip()[1:].split(',')[1:]

data = np.loadtxt('vilar.csv', delimiter=',', skiprows=1)
ns = data[:, 0]

for i in [3, 1, 2, 4, 0]:
    name = header[i]
    times = data[:, i+1]
    res = stats.linregress(ns, times)
    plt.plot(ns, times, 'o', color=f'C{i}', label=name)
    plt.plot(ns, res.intercept + res.slope * ns,
             label=f'{res.intercept:.3f} + {res.slope:.3f} * N', color=f'C{i}')

plt.xlabel('Number of trajectories')
plt.ylabel('Time [s]')
plt.title('Vilar oscillator')
plt.grid(True)
plt.legend()
plt.savefig('vilar.png', bbox_inches='tight')
# plt.show()
