import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import matplotlib.gridspec as gridspec
import h5py

# Set font to monospace since animation looks better this way
plt.rcParams.update({"font.family": "monospace"})


def make_plot(axis, data_file_name):
    # We want to plot the generated data
    hf = h5py.File("./data/" + data_file_name, "r")
    # Names of the planets and keys to access the different datasets
    names = list(hf.keys())
    names = [x for x in hf.keys() if x!="Time" if x!="Pluto"]
    # separate time data from spatial data
    time = np.array(hf.get("Time"), dtype=float)

    # Get the spatial data of each stellar object and save it into a list
    planets = []

    for name in names:
        planets.append(np.array(hf.get(name), dtype=float).transpose())
    planets = np.array(planets)

    # Setting title and other cosmetic parameters for plot
    if data_file_name[-7:-5]=="ee":
        axis.set_title("Explicit Euler")
    elif data_file_name[-7:-5]=="vv":
        axis.set_title("Velocity Verlet")
    elif data_file_name[-7:-5]=="em":
        axis.set_title("Explicit Midpoint")
    else:
        axis.set_title("ODE Int")
    """
    # Change major ticks to show every 20.
    axis.xaxis.set_major_locator(MultipleLocator(20))
    axis.yaxis.set_major_locator(MultipleLocator(20))

    # Change minor ticks to show every 5. (20/4 = 5)
    axis.xaxis.set_minor_locator(AutoMinorLocator(4))
    axis.yaxis.set_minor_locator(AutoMinorLocator(4))

    # Turn grid on for both major and minor ticks and style minor slightly
    # differently.
    axis.grid(which='major', color='#CCCCCC', linestyle='--')
    axis.grid(which='minor', color='#CCCCCC', linestyle=':')
    """
    # Set the initial positions of all the planets from which the animation starts
    lines2d = [
        axis.plot(  planet[0,:-1], 
                    planet[1,:-1], 
                    planet[2,:-1], 
                    linestyle="dashdot", 
                    label=name)[0] 
        for planet, name in zip(planets,names)]

    return lines2d, planets, time


if __name__ == "__main__":
    # Set up environment for plot
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 4)
    gs.update(wspace=0.5)
    ax1 = fig.add_subplot(2,2,1, projection='3d')
    ax2 = fig.add_subplot(2,2,2, projection='3d')
    ax3 = fig.add_subplot(2,2,3, projection='3d')
    ax4 = fig.add_subplot(2,2,4, projection='3d')

    # Get data for plot. Note that the subscript stands for the different 
    # integration schemes:
    #   ee: Explicit Euler
    #   em: Explicit Midpoint
    #   vv: Velocity Verlet
    line_ee, planets_ee, time_ee = make_plot(ax1, "naive_ee.hdf5")
    line_em, planets_em, time_em = make_plot(ax2, "naive_em.hdf5")
    line_vv, planets_vv, time_vv = make_plot(ax3, "naive_vv.hdf5")
    line_od, planets_od, time_od = make_plot(ax4, "ode_int.hdf5")

    # A line in this context is the data that we hand pythons animation function later
    lines = [line_ee, line_em, line_vv, line_od]
    spatial_data_container_list = [planets_ee, planets_em, planets_vv, planets_od]
    time_container = [time_ee, time_em, time_vv, time_od]

    def update_curves(num):
        num *=500
        for line, planets, time in zip(lines, spatial_data_container_list, time_container):
            if num > len(planets[0][0]):
                # Freeze animation at last frame 
                num = len(planets[0][0])-1

            for line, planet in zip(line, planets):
                line.set_data(planet[0, :num], planet[1, :num])
                line.set_3d_properties(planet[2, :num])
        fig.suptitle('September 5th, 1994 at 00:00h\n + {:3}y {:3}d'.format(int(time[num]/365),int(time[num]%365)), y=0.95)
        return lines,

    ani = matplotlib.animation.FuncAnimation(fig, update_curves, frames=3000, interval=60, blit=False, save_count=200)
    for axis in [ax1, ax2, ax3, ax4]:
        axis.view_init(elev=20, azim=100)
        axis.dist=10
        axis.set_xlim(axis.get_xlim())
        axis.set_ylim(axis.get_ylim())
        axis.set_zlim(axis.get_zlim())
        axis.set_xlabel("[A.U.]")
        axis.set_ylabel("[A.U.]")
        axis.set_zlabel("[A.U.]")
    ax1.legend(loc=(0.15,-1.4), ncol=6)
    plt.show()
