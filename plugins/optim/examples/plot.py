import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math


with open("output.txt") as f:
    reader = f.read()
    reader = reader.replace('inf', "math.inf")
    LB = eval(reader.split("\n")[0])
    UB = eval(reader.split("\n")[1])
    UBx = []
    UBy = []
    fig = plt.figure()
    ax1 = fig.add_subplot(111, aspect='equal')
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
    # plt.plot(LBx, LBy, 'r.')
    plt.plot(UBx, UBy, 'b.')
    plt.legend(['UB', 'LB'], loc='upper right')
    plt.show()
