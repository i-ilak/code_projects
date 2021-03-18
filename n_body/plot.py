import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

import matplotlib as mpl
plt.rcParams.update({"font.family": "monospace"})

# We want to plot the generated data form data.txt.

file = np.array(open("data.txt").read().split())
N = int(file[0])    # the first line only contains the number of stellar objects that we simulated.
names = file[1].split(",")
file = np.array(file[2:], dtype=float)

# We then reshape the file to get into the form we are used. 
file = file.reshape((int(len(file)/(6*N+1)), (6*N+1)))

time = file[:,-1]

planets = []

for k in range(N):
    planet = []
    x = file[:, 3*k]
    y = file[:, 3*k+1]
    z = file[:, 3*k+2]
    planet.append(x)
    planet.append(y)
    planet.append(z)
    planets.append(planet)
planets = np.array(planets)

fig = plt.figure()
ax = plt.axes(projection='3d')
title = ax.set_title("September 5th, 1994 at 00:00h\n + 0 d")

# Change major ticks to show every 20.
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_major_locator(MultipleLocator(20))

# Change minor ticks to show every 5. (20/4 = 5)
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))

# Turn grid on for both major and minor ticks and style minor slightly
# differently.
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')

def update_curves(num):
    num *=500
    if num > len(planets[0][0]):
        # Freeze animation at last frame 
        num = len(planets[0][0])-1

    for line, planet in zip(lines2d, planets):
        line.set_data(planet[0, :num], planet[1, :num])
        line.set_3d_properties(planet[2, :num])
    ax.view_init(elev=20, azim=100)
    ax.dist=10
    title.set_text('September 5th, 1994 at 00:00h\n + {:3}y {:3}d'.format(int(time[num]/365),int(time[num]%365)))
    return title, lines2d

lines2d = [ax.plot(planet[0,:-1], planet[1,:-1], planet[2,:-1], linestyle="dashdot", label=name)[0] for planet,name in zip(planets,names)]

ani = matplotlib.animation.FuncAnimation(fig, update_curves, frames=3000, interval=60, blit=False, save_count=200)
ax.set_xlim(ax.get_xlim())
ax.set_ylim(ax.get_ylim())
ax.set_zlim(ax.get_zlim())
ax.set_xlabel("[A.U.]")
ax.set_ylabel("[A.U.]")
ax.set_zlabel("[A.U.]")
plt.legend(loc="best")
#plt.tight_layout()
#ani.save("try_animation.mp4", fps=60)
plt.show()