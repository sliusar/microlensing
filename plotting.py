# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import matplotlib.pyplot as plt
import numpy as np
import re
import matplotlib
# %matplotlib notebook

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn

def density_scatter(x, y, bins=50, xlim=None, ylim=None, filename=None):
    fig = plt.figure(figsize=(10, 10))
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.clf()
    plt.imshow(heatmap.T, extent=extent, origin='lower')
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if filename is not None:
        plt.savefig(filename)
    plt.show()
    
import matplotlib.animation as animation
# -





rays_y = np.loadtxt("rays_y_0.00.dat")

m = np.logical_and(np.abs(rays_y[:,0]) <= 20, np.abs(rays_y[:,1]) <= 20)
x = rays_y[:,0][m]
y = rays_y[:,1][m]
density_scatter(x, y, bins=[200, 200])

m = np.logical_and(np.abs(rays_y[:,0]) <= 15, np.abs(rays_y[:,1]) <= 15)
x = rays_y[:,0][m]
y = rays_y[:,1][m]
density_scatter(x, y, bins=[200, 200])

m = np.logical_and(np.abs(rays_y[:,0]) <= 30, np.abs(rays_y[:,1]) <= 30)
x = rays_y[:,0][m]
y = rays_y[:,1][m]
density_scatter(x, y, bins=[500, 500])


def get_image_data(filename, gamma=0.6):
    image = np.loadtxt(filename).reshape([500,500])
    extent = None
    with open(filename, 'r') as f:
        header1 = f.readline().strip()
        header2 = f.readline().strip()

        _x = re.split('\(|\,|\)', header1)
        _y = re.split('\(|\,|\)', header2)
        extent = [float(_x[1]),float(_x[2]),float(_y[1]),float(_y[2])]

    img = np.zeros_like(image)
    img[image > 0] = image[image > 0]**gamma
    if img[img == np.max(img)].shape[0] < 5:
        img[img == np.max(img)] = np.mean(img) # removing center-of-mass pixel with extreame amplification
    return img, extent
    


# +
img, extent = get_image_data("images/image_0.00.dat", gamma=1)
fig, ax = plt.subplots(figsize=(10,10))
pos = ax.imshow(img, interpolation='none', extent=extent, origin='lower')
fig.colorbar(pos, ax=ax)
plt.show()

img, extent = get_image_data("images/image_0.90.dat", gamma=1)
fig, ax = plt.subplots(figsize=(10,10))
pos = ax.imshow(img, interpolation='none', extent=extent, origin='lower')
fig.colorbar(pos, ax=ax)
plt.show()
# -



a=1

# +
fig, ax = plt.subplots(figsize=(10,10))

ims = []

for i in np.arange(0, 1, 0.1):
    filename = "images/image_%.2f.dat" % i
    img, extent = get_image_data(filename, gamma=1)
    title = plt.text(0.5,1.01, "t=%.2f" % i, ha="center",va="bottom", transform=ax.transAxes, fontsize="large")
    text = ax.text('','','')
    im = ax.imshow(img, extent=extent, origin='lower')
    ims.append([text, im, title])

ani = animation.ArtistAnimation(fig, ims, interval=200, blit=False, repeat=False)
ani.save("images/moving_stars.mp4")
plt.show()
# -


































