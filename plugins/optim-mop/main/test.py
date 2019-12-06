from __future__ import print_function
"""
Do a mouseclick somewhere, move the mouse to some destination, release
the button.  This class gives click- and release-events and also draws
a line or a box from the click-point to the actual mouseposition
(within the same axes) until the button is released.  Within the
method 'self.ignore()' it is checked whether the button from eventpress
and eventrelease are the same.

"""
import numpy as np

import matplotlib.backends.backend_pdf
import time

import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib import collections  as mc
import math
import os
import random
from stat import *
import multiprocessing
from tkinter import *
import time
from itertools import islice
from matplotlib.widgets import Button
# window=Tk()
filedate = None

global filedate
# filedate = None
q = multiprocessing.Queue()

#Create and start the simulation process
simulate=multiprocessing.Process(None,simulation,args=(q,))
simulate.start()
  


#Create the base plot
plot()


#Call a function to update the plot when there is new data
updateplot(q)

print('Done')
    
def line_select_callback(eclick, erelease):
    f= open("selectedSpace.txt","w+")
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    f.write("%3.2f,%3.2f,%3.2f,%3.2f" % (x1, y1, x2, y2))
    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
    print(" The button you used were: %s %s" % (eclick.button, erelease.button))


def toggle_selector(event):
    print(' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print(' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)
        
def _yes(event):
    f= open("selectedSpace.txt","w+")
    f.write("")

print("\n      click  -->  release")

def plot():
    global line,ax1,canvas
    fig = plt.figure()
    plt.subplots_adjust(bottom=0.2)
    axcut = plt.axes([0.7, 0.05, 0.1, 0.075])
    bcut = Button(axcut, 'Terminar',color='red', hovercolor='green')
    ax1 = fig.add_subplot(1,1,1)
    bcut.on_clicked(_yes)
    # plt.ion()
    # canvas = FigureCanvasTkAgg(fig, master=window)
    # canvas.show()
    # canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    # canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    line, = ax1.plot([], [])
    # drawtype is 'box' or 'line' or 'none'
    toggle_selector.RS = RectangleSelector(ax1, line_select_callback,
                                           drawtype='box', useblit=True,
                                           button=[1, 3],  # don't use middle button
                                           minspanx=5, minspany=5,
                                           spancoords='pixels',
                                           interactive=True)
    plt.connect('key_press_event', toggle_selector)

def updateplot(q):
    try:       #Try to check if there is data in the queue
        result=q.get_nowait()

        if result !='Q':
            # print(result)
            LB, UB, LB2 = result
            UBx = []
            UBy = []
            LBx = []
            LBy = []

            ax1.clear()

            for ub in UB:
                UBx.append(ub[0])
                UBy.append(ub[1])

            # add line upperbound
            lines = []
            #p = [(-100, 100), (100, -100)]

            if len(LB)>0:
                lines.append(LB)
                lc = mc.LineCollection(lines, colors='blue', linewidths=0.5)
                ax1.add_collection(lc)

            # add lines function
            lines = []
            # p = [(100,100),(20,60),(-100,-100)]
            if len(LB2)>0:
                for i in range(1,len(LB2)):
                    arr1 = []
                    arr1.append(LB2[i-1])
                    arr1.append(LB2[i])
                    lines.append(arr1)
                lc = mc.LineCollection(lines, colors='green', linewidths=0.5)
                ax1.add_collection(lc)
            ax1.plot()
            plt.plot(UBx, UBy, 'r-', markersize=3)
            # plt.plot(LBx, LBy, '-b', lw=0.5)
            #ax1.set_xlim([0.0,0.35])
            #ax1.set_ylim([0.7,1.0])
            #plt.gca().set_aspect('equal', adjustable='box')
            plt.axis('scaled')
            plt.pause(1)
            
            #if os.path.exists('result.png'):
            #    plt.savefig('result_{}.png'.format(int(time.time())))
            #else:
            #    plt.savefig('result.png')
            
            updateplot(q)
             # print(result)
                 #here get crazy with the plotting, you have access to all the global variables that you defined in the plot function, and have the data that the simulation sent.
             # line.set_ydata([1,result,10])
             # ax1.draw_artist(line)
             # plt.draw()
             # plt.pause(0.1)
             # updateplot(q)
             # window.after(500,updateplot,q)
        else:
             print('done')
    except:
        # print("empty")
        plt.pause(1)
        updateplot(q)
        # window.after(500,updateplot,q)

def simulation(q):
    global filedate
    # print(filedate)
    while True:
        st = os.stat('output2.txt')
        # print(st[ST_MTIME])
        if(st[ST_MTIME] == filedate):
            time.sleep(1)
            pass
        else:
            filedate = st[ST_MTIME]
            f = open("output2.txt")
            reader = f.read()
            reader = reader.replace('inf', "math.inf")
            reader = reader.replace('nan', "0")
            UB = eval(reader.split("\n")[0])
            LB = eval(reader.split("\n")[1])
            LB2 = [] #eval(reader.split("\n")[2])

            # LB2 = eval(reader.split("\n")[2])
            q.put((LB, UB, LB2))

    # iterations = range(100)
    # for i in iterations:
    #     if not i % 10:
    #         time.sleep(1)
    #             #here send any data you want to send to the other process, can be any pickable object
    #         q.put(random.randint(1,10))
    # q.put('Q')


#fig, current_ax = plt.subplots()                 # make a new plotting range
#N = 100000                                       # If N is large one can see
#x = np.linspace(0.0, 10.0, N)                    # improvement by use blitting!

#plt.plot(x, +np.sin(.2*np.pi*x), lw=3.5, c='b', alpha=.7)  # plot something
#plt.plot(x, +np.cos(.2*np.pi*x), lw=3.5, c='r', alpha=.5)
#plt.plot(x, -np.sin(.2*np.pi*x), lw=3.5, c='g', alpha=.3)



#axes = plt.gca()
#axes.set_xlim([0,20])
#axes.set_ylim([0,20])
