import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

U = pa.cx_cube()
U.load("./data/system2_U.bin")
U = np.array(U).T
print(U.shape)
# plt.imshow(U[0,:,:].imag)



# Plot conjugate of U
# plt.imshow(U[0,:,:].conj().imag)


P = np.sqrt(np.multiply(U[:,:,0], U[:,:,0].conj()))
plt.imshow(P.real)
plt.colorbar()
plt.show()

P = np.sqrt(np.multiply(U[:,:,19], U[:,:,19].conj()))
plt.imshow(P.real)
plt.colorbar()
plt.show()

P = np.sqrt(np.multiply(U[:,:,29], U[:,:,29].conj()))
plt.imshow(P.real)
plt.colorbar()
plt.show()


P = np.sqrt(np.multiply(U[:,:,39], U[:,:,39].conj()))
plt.imshow(P.real)
plt.colorbar()
plt.show()

P = np.sqrt(np.multiply(U[:,:,49], U[:,:,49].conj()))
plt.imshow(P.real)
plt.colorbar()
plt.show()

P = np.sqrt(np.multiply(U[:,:,59], U[:,:,59].conj()))
plt.imshow(P.real)
plt.colorbar()
plt.show()

P = np.sqrt(np.multiply(U[:,:,69], U[:,:,69].conj()))
plt.imshow(P.real)
plt.colorbar()
plt.show()

P = np.sqrt(np.multiply(U[:,:,79], U[:,:,79].conj()))
plt.imshow(P.real)
plt.colorbar()
plt.show()

V = pa.mat()
V.load("./data/system2_V.bin")
V = np.array(V).T

plt.imshow(V)
plt.colorbar()
plt.show()



