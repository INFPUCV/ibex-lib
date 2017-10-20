import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import math


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
        for p in [patches.Rectangle(
                    (lb['pts'][0], lb['pts'][1]),
                    lb['diam_x'], lb['diam_y'],
                    fill=False,
                    edgecolor="red",
                    linestyle='dashed'
                    ) for lb in LB]:
            ax1.add_patch(p)
        for ub in UB:
            UBx.append(ub[0])
            UBy.append(ub[1])
        ax1.plot()
        plt.legend(['UB', 'LB'], loc='upper right')
        ax1.plot(UBx, UBy, 'b.')
        # plt.plot(LBx, LBy, 'r.')
        # plt.show()
    except SyntaxError:
        print('holi')

ani = animation.FuncAnimation(fig, animate)
plt.show()
