import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import ticker
import pyarma as pa

# Load U
U = pa.cx_cube()
U.load("./data/system2_U.bin")
U = np.array(U).T

# Load V
V = pa.mat()
V.load("./data/system2_V.bin")
V = np.array(V).T

# Remove parts of V that are not the barrier for transparency
V[V < 100] = np.nan

# Fontsize settings
label_fontsize = 23
ticks_fontsize = 20
tick_font_size = 17


# PLOT PROBABILITY WAVE FUNCTION

# Plot probability wave function at t = 0.000
fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,0], U[:,:,0].conj()))
plt.imshow(P.real, extent=[0,1,0,1])
plt.imshow(V,extent=[0,1,0,1], cmap = ListedColormap(['black']))
plt.xlabel(r'$x$', fontsize=label_fontsize+2)
plt.ylabel(r'$y$', fontsize=label_fontsize+2)
plt.xticks(fontsize=ticks_fontsize, rotation=40)
plt.yticks(fontsize=ticks_fontsize, rotation=40)
plt.title("t = 0.000", fontsize=label_fontsize)
ax = plt.gca()
ax.xaxis.set_ticks([0, 0.5, 1])
ax.yaxis.set_ticks([0, 0.5, 1])
plt.tight_layout()
plt.savefig(f"./figs/system2_U_000.pdf")

# Plot probability wave function at t = 0.001
fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,39], U[:,:,39].conj()))
plt.imshow(P.real, extent=[0,1,0,1])
plt.imshow(V,extent=[0,1,0,1], cmap = ListedColormap(['black']))
plt.xlabel(r'$x$', fontsize=label_fontsize+2)
plt.xticks(fontsize=ticks_fontsize, rotation=40)
plt.yticks(fontsize=ticks_fontsize, rotation=40)
plt.title("t = 0.001", fontsize=label_fontsize)
ax = plt.gca()
ax.xaxis.set_ticks([0, 0.5, 1])
ax.yaxis.set_ticks([0, 0.5, 1])
plt.tight_layout()
plt.savefig(f"./figs/system2_U_001.pdf")

# Plot probability wave function at t = 0.002
fig = plt.figure(figsize=(6, 4.5))
P = np.sqrt(np.multiply(U[:,:,79], U[:,:,79].conj()))
im = plt.imshow(P.real, extent=[0,1,0,1])
plt.imshow(V,extent=[0,1,0,1], cmap = ListedColormap(['black']))
plt.xlabel(r'$x$', fontsize=label_fontsize+2)
plt.xticks(fontsize=ticks_fontsize, rotation=40)
plt.yticks(fontsize=ticks_fontsize, rotation=40)
plt.title("t = 0.002", fontsize=label_fontsize)
ax = plt.gca()
ax.xaxis.set_ticks([0, 0.5, 1])
ax.yaxis.set_ticks([0, 0.5, 1])
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=tick_font_size)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.tight_layout()
plt.savefig(f"./figs/system2_U_002.pdf")
cbar.remove()

# Reset matplotlib
plt.rcParams.update(plt.rcParamsDefault)


# REAL PART #

# Plot real part of wave function at t = 0.000
fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,0].real
plt.imshow(P, extent=[0,1,0,1])
plt.imshow(V,extent=[0,1,0,1], cmap = ListedColormap(['black']))
plt.xlabel(r'$x$', fontsize=label_fontsize+2)
plt.ylabel(r'$y$', fontsize=label_fontsize+2)
plt.xticks(fontsize=ticks_fontsize, rotation=40)
plt.yticks(fontsize=ticks_fontsize, rotation=40)
plt.title("t = 0.000 (real)", fontsize=label_fontsize)
ax = plt.gca()
ax.xaxis.set_ticks([0, 0.5, 1])
ax.yaxis.set_ticks([0, 0.5, 1])
plt.tight_layout()
plt.savefig(f"./figs/system2_U_000_r.pdf")

# Plot real part of wave function at t = 0.001
fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,39].real
plt.imshow(P, extent=[0,1,0,1])
plt.imshow(V,extent=[0,1,0,1], cmap = ListedColormap(['black']))
plt.xlabel(r'$x$', fontsize=label_fontsize+2)
plt.xticks(fontsize=ticks_fontsize, rotation=40)
plt.yticks(fontsize=ticks_fontsize, rotation=40)
plt.title("t = 0.001 (real)", fontsize=label_fontsize)
ax = plt.gca()
ax.xaxis.set_ticks([0, 0.5, 1])
ax.yaxis.set_ticks([0, 0.5, 1])
plt.tight_layout()
plt.savefig(f"./figs/system2_U_001_r.pdf")

# Plot real part of wave function at t = 0.002
fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,79].real
im = plt.imshow(P, extent=[0,1,0,1])
plt.imshow(V,extent=[0,1,0,1], cmap = ListedColormap(['black']))
plt.xlabel(r'$x$', fontsize=label_fontsize+2)
plt.xticks(fontsize=ticks_fontsize, rotation=40)
plt.yticks(fontsize=ticks_fontsize, rotation=40)
plt.title("t = 0.002 (real)", fontsize=label_fontsize)
ax = plt.gca()
ax.xaxis.set_ticks([0, 0.5, 1])
ax.yaxis.set_ticks([0, 0.5, 1])
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=tick_font_size)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.tight_layout()
plt.savefig(f"./figs/system2_U_002_r.pdf")
cbar.remove()

# Reset matplotlib
plt.rcParams.update(plt.rcParamsDefault)


# IMAGINARY PART #

# Plot imaginary part of wave function at t = 0.000
fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,0].imag
plt.imshow(P, extent=[0,1,0,1])
plt.imshow(V,extent=[0,1,0,1], cmap = ListedColormap(['black']))
plt.xlabel(r'$x$', fontsize=label_fontsize+2)
plt.ylabel(r'$y$', fontsize=label_fontsize+2)
plt.xticks(fontsize=ticks_fontsize, rotation=40)
plt.yticks(fontsize=ticks_fontsize, rotation=40)
plt.title("t = 0.000 (imag)", fontsize=label_fontsize)
ax = plt.gca()
ax.xaxis.set_ticks([0, 0.5, 1])
ax.yaxis.set_ticks([0, 0.5, 1])
plt.tight_layout()
plt.savefig(f"./figs/system2_U_000_i.pdf")


# Plot imaginary part of wave function at t = 0.001
fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,39].imag
plt.imshow(P, extent=[0,1,0,1])
plt.imshow(V,extent=[0,1,0,1], cmap = ListedColormap(['black']))
plt.xlabel(r'$x$', fontsize=label_fontsize+2)
plt.xticks(fontsize=ticks_fontsize, rotation=40)
plt.yticks(fontsize=ticks_fontsize, rotation=40)
plt.title("t = 0.001 (imag)", fontsize=label_fontsize)
ax = plt.gca()
ax.xaxis.set_ticks([0, 0.5, 1])
ax.yaxis.set_ticks([0, 0.5, 1])
plt.tight_layout()
plt.savefig(f"./figs/system2_U_001_i.pdf")

# Plot imaginary part of wave function at t = 0.002
fig = plt.figure(figsize=(6, 4.5))
P = U[:,:,79].imag
im = plt.imshow(P, extent=[0,1,0,1])
plt.imshow(V,extent=[0,1,0,1], cmap = ListedColormap(['black']))
plt.xlabel(r'$x$', fontsize=label_fontsize+2)
plt.xticks(fontsize=ticks_fontsize, rotation=40)
plt.yticks(fontsize=ticks_fontsize, rotation=40)
plt.title("t = 0.002 (imag)", fontsize=label_fontsize)
ax = plt.gca()
ax.xaxis.set_ticks([0, 0.5, 1])
ax.yaxis.set_ticks([0, 0.5, 1])
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=tick_font_size)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.tight_layout()
plt.savefig(f"./figs/system2_U_002_i.pdf")
cbar.remove()


