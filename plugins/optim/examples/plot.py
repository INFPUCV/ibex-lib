import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import math

from itertools import islice

fig = plt.figure()
ax1 = fig.add_subplot(111, aspect='equal')

def animate(i):
    try:
        f = open("output.txt")
        reader = f.read()
        reader = reader.replace('inf', "math.inf")
        LB = eval(reader.split("\n")[0])
        UB = eval(reader.split("\n")[1])
        UBx = []
        UBy = []
        ax1.clear()

        ax1.add_patch(patches.Rectangle(
            (LB[0]['pts'][0], LB[0]['pts'][1]),
            LB[0]['diam_x'], LB[0]['diam_y'],
            fill=False,
            edgecolor='red',
            linestyle='solid',
            lw=0.1
            ))

        for p in [patches.Rectangle(
                    (lb['pts'][0], lb['pts'][1]),
                    lb['diam_x'], lb['diam_y'],
                    fill=False,
                    edgecolor='black',
                    linestyle='solid',
                    lw=0.1
                    ) for lb in islice(LB, 1, len(LB))]:
            ax1.add_patch(p)
        for ub in UB:
            UBx.append(ub[0])
            UBy.append(ub[1])
        for lb in LB:
        	if (lb['pA'][0] < lb['pB'][0]) and (lb['pA'][1] > lb['pB'][1]):
        		line = plt.Line2D((lb['pA'][0], lb['pB'][0]),(lb['pA'][1], lb['pB'][1]),lw=0.5,markeredgecolor='black')
        		ax1.add_line(line)
        ax1.plot()
        plt.legend(['UB', 'LB'], loc='upper right')
        # ax1.plot(UBx, UBy, 'b.')
        plt.plot(UBx, UBy, 'r.', markersize=1)

        #plt.show()
    except SyntaxError:
        print('holi')

ani = animation.FuncAnimation(fig, animate)
plt.show()
