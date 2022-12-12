#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pyarma as pa

# Load data
U = pa.cx_cube()
U.load("./data/system2_U.bin")
U = np.array(U).T

# Calculate probability and change to real matrix
U = np.multiply(U, U.conj())
U = U.real

# Create figure adn remove white space
fig = plt.figure()
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)


# Plot the first frame
img = ax.imshow(U[:,:,0], cmap=plt.get_cmap("viridis"))

# Remove axis ticks
ax.set_xticks([])
ax.set_yticks([])

# Animation function to update the plot
def animation(i):
    img.set_data(U[:,:,i])

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, 80, 2), repeat=False, blit=0)

# Save the animation as a GIF
anim.save("./figs/U.gif", writer="imagemagick", fps=24)

