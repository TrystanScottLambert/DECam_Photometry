"""
Toy model of a filament to see if we can recover the shape that we are seeing due to the 
offset in narrowband filter.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import norm


number = 1000
elevation = 90 # degrees
azim = -20 # degrees
r = np.random.random(number) * 150
theta = np.random.random(number) * np.pi * 2

y = np.random.random(number) * 10
x = np.sqrt(r) * np.cos(theta)
z = np.sqrt(r) * np.sin(theta)


fig = plt.figure()
ax = fig.add_subplot(121,projection='3d')
ax.scatter(x, y, z)
ax.view_init(elev = elevation, azim=azim)

ax_proj = fig.add_subplot(122)

points = np.vstack([x, y, z]).T

azim = ax.azim*np.pi/180
elev = ax.elev*np.pi/180
elev *= 1.2          # this seems to improve the outcome

a_vec = np.array([np.cos(azim),np.sin(azim),0])
normal = np.cos(elev)*a_vec + np.array([0,0,np.sin(elev)])

z_vec = np.array([0,0,1])
y_comp = z_vec - (z_vec@normal)*normal
y_comp = y_comp/norm(y_comp)
x_comp = np.cross(y_comp,normal)

proj_mat = np.vstack([x_comp,y_comp]) # build projection matrix
proj_mat = -proj_mat                  # account for flipped axes
points_2D = points @ proj_mat.T       # apply projection

ax_proj.scatter(*points_2D.T)
#plt.gca().set_aspect('equal', adjustable='box')
plt.show()

