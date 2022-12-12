#!/usr/bin/env python3
import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation, rc
from IPython.display import HTML, Image

#
# Let's generate a dummy time series for a function z(x,y,t)
#

# Set up a 2D xy grid
# h = 0.005
# x_points = np.arange(0, 1+h, h)
# y_points = np.arange(0, 1+h, h)
# x, y = np.meshgrid(x_points, y_points, sparse=True)

# # Array of time points
# dt = 0.005
# t_points = np.arange(0, 1+dt, dt)

# # A function for a Gaussian that is travelling
# # in the x direction and broadening as time passes
# def z(x,y,t):
#     v = 0.5
#     x_c = 0.2
#     sigma_x = 0.025 + 0.15 * t
#     return 1. / (2 * np.pi * np.sqrt(sigma_x)) * np.exp(-0.5 * (x - x_c - v * t)**2 / sigma_x**2)

U = pa.cx_cube()
U.load("./data/system2_U.bin")
U = np.array(U).T

U = np.multiply(U, U.conj())
U = U.real

t_points = 80




# Fill z_data_list with f(x,y,t)
z_data_list = []
for t in range(t_points):
    z_data = U[:,:,t]
    z_data_list.append(z_data)


#
# Now the list z_data_list contains a series of "frames" of z(x,y,t),
# where each frame can be plotted as a 2D image using imshow. Let's
# animate it!
#

# Some settings
fontsize = 12
# t_min = t_points[0]
# x_min, x_max = x_points[0], x_points[-1]
# y_min, y_max = y_points[0], y_points[-1]

# Create figure
fig = plt.figure()
ax = plt.gca()

# Create a colour scale normalization according to the max z value in the first frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

# Plot the first frame
img = ax.imshow(U[:,:,0], cmap=plt.get_cmap("viridis"), norm=norm)
##Remove axis ticks
ax.set_xticks([])
ax.set_yticks([])



# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}", color="white",
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
    img.set_norm(norm)

    # Update z data
    img.set_data(z_data_list[i])

    # Update the time label
    current_time = 0 + i * 0.01
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=False, blit=0)

# Run the animation!
plt.show()
anim.save("./figs/U.gif", writer="imagemagick", fps=24)

# # Save the animation
# anim.save('./animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)
