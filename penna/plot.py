import numpy as np
import matplotlib.pyplot as plt

file = np.array(open("data.txt").read().split(), dtype=int)
file = file.reshape((int(len(file)/2), 2))

plt.figure()
plt.loglog(file[:, 0], file[:, 1], "r")
plt.xlabel(r"Generation")
plt.ylabel(r"$n$")
plt.tight_layout()
plt.savefig("pop.pdf")
