import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

U = pa.cx_cube()
U.load("./data/system2_U.bin")
U = np.array(U).T
print(U.shape)
#
#
# TODO: NORMALIZE AXIS
#       ADD PROBABILITY PLOT
#       CLEAN UP CODE
#

# Fontsize
label_fontsize = 16
ticks_fontsize = 16
legend_fontsize = 16


# Real part
fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,0], U[:,:,0].conj()))
plt.imshow(P.real)
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_000_r.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,39], U[:,:,39].conj()))
plt.imshow(P.real)
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_001_r.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,79], U[:,:,79].conj()))
plt.imshow(P.real)
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_002_r.pdf")

# Imaginary part
fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,0], U[:,:,0].conj()))
plt.imshow(P.imag)
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_000_i.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,39], U[:,:,39].conj()))
plt.imshow(P.imag)
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize) 


plt.savefig(f"./figs/system2_U_001_i.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,79], U[:,:,79].conj()))
plt.imshow(P.imag)
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_002_i.pdf")


# Plot potential
V = pa.mat()
V.load("./data/system2_V.bin")
V = np.array(V).T

plt.imshow(V)
plt.colorbar()
plt.savefig("./figs/system2_V.pdf")



