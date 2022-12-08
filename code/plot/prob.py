import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

P = pa.mat()
P.load("./data/system1_prob.bin")
P = np.array(P)

n = len(P)

plt.plot(range(n), P, 'b-')
plt.xlabel('n')
plt.ylabel('P(n)')
plt.show()