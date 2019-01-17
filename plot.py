from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.animation as animation

figSalida = plt.figure()
currentAxis = plt.gca()

def setDim():
    with open("outDim.txt") as f:
        content = f.readlines()
        ly = content[0].strip("[").strip("]\n").strip(" ").split(",")
        lx = content[1].strip("[").strip("]\n").strip(" ").split(",")

    currentAxis.set_xlim([0, float(lx[0])])
    currentAxis.set_ylim([0, float(ly[0])])

setDim()
rects = []

def anim(i):
    #setDim()
    global rects
    for rect in rects:
        rect.remove()
    rects = []

    with open("cajasCurrent.txt") as f:
        content = f.readlines()

        for index in range(int(len(content)/2)):
            index = index*2
            ly = content[index].strip("[").strip("]\n").strip(" ").split(",")
            index += 1
            lx = content[index].strip("[").strip("]\n").strip(" ").split(",")

            pos = (float(lx[0]), float(ly[0]))
            dimX = float(lx[1]) - float(lx[0])
            dimY = float(ly[1]) - float(ly[0])

            rects.append(Rectangle(pos, dimX, dimY, 
                alpha=0.3, facecolor="green", edgecolor="black"))
            currentAxis.add_patch(rects[-1])
    
    with open("cajasDescartadas.txt") as f:
        content = f.readlines()

        for index in range(int(len(content)/2)):
            index = index*2
            ly = content[index].strip("[").strip("]\n").strip(" ").split(",")
            index += 1
            lx = content[index].strip("[").strip("]\n").strip(" ").split(",")

            pos = (float(lx[0]), float(ly[0]))
            dimX = float(lx[1]) - float(lx[0])
            dimY = float(ly[1]) - float(ly[0])

            rects.append(Rectangle(pos, dimX, dimY, 
                alpha=0.3, facecolor="red", edgecolor="black"))
            currentAxis.add_patch(rects[-1])
    

ani = animation.FuncAnimation(figSalida, anim, interval=1000)
plt.show()
