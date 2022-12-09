import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

U = pa.mat()
U.load("./data/system1_U.bin")
U = np.array(U)

# plt.imshow(U, cmap='hot', interpolation='nearest')
# plt.show()