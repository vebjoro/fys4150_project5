import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa



# Fontsize
label_fontsize = 16
ticks_fontsize = 16


# 2 SLIT
U = pa.cx_cube()
U.load("./data/system2_U.bin")
U = np.array(U).T

U = np.multiply(U, U.conj())

U_slice = U[:,:,-1]
U_slice = U_slice.real
print(U_slice.shape)
indx = int(np.round(0.8*len(U_slice)))

tot = 0
for i in range(len(U_slice)):
    tot += U_slice[i,indx]

fig = plt.figure(figsize=(6, 4.5))
plt.plot(range(len(U_slice)), U_slice[:,indx]/tot, color="#8a1629", alpha=0.8, linewidth=2.0)
plt.xlabel(r"$y$", fontsize=label_fontsize)
plt.ylabel(r"$p(y \mid x=0.8 ; t=0.002)$", fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.tight_layout()
plt.savefig("./figs/2_slit.pdf")



# 1 SLIT
U = pa.cx_cube()
U.load("./data/system2.1_U.bin")
U = np.array(U).T

U = np.multiply(U, U.conj())

U_slice = U[:,:,-1]
U_slice = U_slice.real
print(U_slice.shape)
indx = int(np.round(0.8*len(U_slice)))

tot = 0
for i in range(len(U_slice)):
    tot += U_slice[i,indx]

fig = plt.figure(figsize=(6, 4.5))
plt.plot(range(len(U_slice)), U_slice[:,indx]/tot, color="#8a1629", alpha=0.8, linewidth=2.0)
plt.xlabel(r"$y$", fontsize=label_fontsize)
plt.ylabel(r"$p(y \mid x=0.8 ; t=0.002)$", fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.tight_layout()
plt.savefig("./figs/1_slit.pdf")



# 3 SLIT
U = pa.cx_cube()
U.load("./data/system2.3_U.bin")
U = np.array(U).T

U = np.multiply(U, U.conj())

U_slice = U[:,:,-1]
U_slice = U_slice.real
print(U_slice.shape)
indx = int(np.round(0.8*len(U_slice)))

tot = 0
for i in range(len(U_slice)):
    tot += U_slice[i,indx]

fig = plt.figure(figsize=(6, 4.5))
plt.plot(range(len(U_slice)), U_slice[:,indx]/tot, color="#8a1629", alpha=0.8, linewidth=2.0)
plt.xlabel(r"$y$", fontsize=label_fontsize)
plt.ylabel(r"$p(y \mid x=0.8 ; t=0.002)$", fontsize=label_fontsize)
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)
plt.grid()
ax = plt.gca()
ax.set_facecolor("#e6e6e6")
plt.tight_layout()
plt.savefig("./figs/3_slit.pdf")