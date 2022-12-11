import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

U = pa.cx_cube()
U.load("./data/system2_U.bin")
U = np.array(U).T
print(U.shape)
#
#
# TODO: ADD COLORBAR (trengs kun på siste del av tidsserien t = 0.002)
#       ADD PROBABILITY PLOT
#       CLEAN UP CODE
#       SHOW WALL
#

# Fontsize
label_fontsize = 16
ticks_fontsize = 16
legend_fontsize = 16



fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,0], U[:,:,0].conj()))
plt.imshow(P.real, extent=[0,1,0,1])
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_000.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,39], U[:,:,39].conj()))
plt.imshow(P.real, extent=[0,1,0,1])

plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)

plt.savefig(f"./figs/system2_U_001.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,79], U[:,:,79].conj()))
plt.imshow(P.real, extent=[0,1,0,1])
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)

plt.savefig(f"./figs/system2_U_002.pdf")


# Real part
fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,0].real
plt.imshow(P, extent=[0,1,0,1])
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_000_r.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,39].real
plt.imshow(P, extent=[0,1,0,1])

plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)

plt.savefig(f"./figs/system2_U_001_r.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,79].real
plt.imshow(P, extent=[0,1,0,1])
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_002_r.pdf")

# Imaginary part
fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,0].imag
plt.imshow(P, extent=[0,1,0,1])
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_000_i.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,39].imag
plt.imshow(P, extent=[0,1,0,1])
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize) 


plt.savefig(f"./figs/system2_U_001_i.pdf")

fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,79].imag
plt.imshow(P, extent=[0,1,0,1])
plt.xlabel('x', fontsize=label_fontsize)
plt.ylabel('y', fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)


plt.savefig(f"./figs/system2_U_002_i.pdf")


# Plot potential
V = pa.mat()
V.load("./data/system2_V.bin")
V = np.array(V).T

plt.imshow(V,extent=[0,1,0,1])
plt.colorbar()
plt.savefig("./figs/system2_V.pdf")



