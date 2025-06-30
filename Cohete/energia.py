import matplotlib.pyplot as plt
import numpy as np

energia = np.loadtxt("Energia.txt", dtype=float, delimiter='\t', skiprows=1).transpose()
# plt.plot(energia[0], 'g-')
# plt.plot(energia[1], 'b-')
# plt.plot(energia[0] + energia[1], 'r-')
plt.plot(energia[2], 'k-')
plt.show()