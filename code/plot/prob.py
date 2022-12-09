import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

P = pa.mat()
P.load("./data/system1_prob.bin")
P = np.array(P).T[0]

n = len(P)
print(P)

plt.plot(range(n), P, 'b-')
plt.xlabel('n')
plt.ylabel('P(n)')
plt.ylim(0,2)
plt.show()