import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation


file = np.array(open("data.txt").read().split(), dtype=float)
N = int(file[0])
file = file[1:]
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
title = ax.set_title("TEST")

def update_curves(num):
    num *=10
    if num > len(planets[0][0]):
        ani.event_source.stop()

    for line, planet in zip(lines2d,planets):
        line.set_data(planet[0, :num], planet[1, :num])
        line.set_3d_properties(planet[2, :num])
    ax.view_init(elev=3, azim=360+(num/100)%(num/100+1))
    title.set_text('{}'.format(time[num]))
    return title, lines2d

lines2d = [ax.plot(planet[0,:1], planet[1,:1], planet[2,:1], linestyle="--")[0] for planet in planets]

ani = matplotlib.animation.FuncAnimation(fig, update_curves, frames=3000, interval=60, blit=False, save_count=200)
ax.set_xlim(ax.get_xlim())
ax.set_ylim(ax.get_ylim())
ax.set_zlim(ax.get_zlim())
ax.axis('off')
plt.tight_layout()
#ani.save("try_animation.mp4", fps=60)
plt.show()