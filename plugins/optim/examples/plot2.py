import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.backend_bases import NavigationToolbar2, Event
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)
import math
import os
import random
from stat import *
import multiprocessing
from tkinter import *
import time
from itertools import islice
import subprocess


class Box:
    def __init__(self, box_string):
        box_string = box_string.replace('inf', "math.inf")
        box_string = box_string.replace('nan', "0")
        box_var = eval(box_string)
        self.createBox(box_var['id'], box_var['pts'][0], box_var['pts'][1],
                       box_var['diam_x'], box_var['diam_y'])

    def __repr__(self):
        return str((self.x, self.y))

    def createBox(self, id_box, x, y,  diam_x, diam_y):
        self.id_box = id_box
        self.x = x
        self.y = y
        self.diam_x = diam_x
        self.diam_y = diam_y


class LB:
    def __init__(self, lb_string):
        lb_string = lb_string.replace('inf', 'math.inf')
        lb_string = lb_string.replace('nan', '0')
        lb_var = eval(lb_string)
        self.key = lb_var['id']
        self.x = lb_var['pts'][0]
        self.y = lb_var['pts'][1]

    def __lt__(self, other):
        if isinstance(other, LB):
            if self.x < other.x:
                return True
            elif self.x == other.x:
                return self.y > other.y
            else:
                return False
        return NotImplemented


class UB:
    def __init__(self, ub_string):
        ub_string = ub_string.replace('inf', 'math.inf')
        ub_string = ub_string.replace('nan', '0')
        ub_var = eval(ub_string)
        self.key = ub_string
        self.x = ub_var['pts'][0]
        self.y = ub_var['pts'][1]


class DataDict:
    box_set = {}
    patch_set = {}
    UB_set = {}
    LB_set = {}
    pending_box = None

    def get_ub_x(self):
        box_list = []
        for x in self.UB_set.values():
            box_list.append(x.x)
        return box_list

    def get_ub_y(self):
        box_list = []
        for x in self.UB_set.values():
            box_list.append(x.y)
        return box_list

    def get_lb(self):
        lb_list_x = []
        lb_list_y = []
        box_list = sorted(self.LB_set.values())
        print(box_list)
        i = 0
        for box in box_list[1:]:
            if box.y < box_list[i].y:
                # print(box.y)
                # print(box_list[i].y)
                box_list.remove(box)
            i = i + 1
        i = 0
        print(box_list)
        # lb_list_x.append(box_list[0].x)
        # lb_list_y.append(box_list[0].y)
        # for box in box_list[1:]:
        #     lb_list_x.append(box.x)
        #     lb_list_y.append(box_list[i].y)
        #     lb_list_x.append(box.x)
        #     lb_list_y.append(box.y)
        #     i = i + 1
        # print(lb_list_x)
        return lb_list_x, lb_list_y

    def append_lb(self, lb):
        self.LB_set[lb.key] = lb

    def append_ub(self, ub):
        self.UB_set[ub.key] = ub

    def remove_ub(self, key):
        key = key.replace('inf', 'math.inf')
        key = key.replace('nan', '0')
        del self.UB_set[key]

    def append_box(self, box):
        self.box_set[box.id_box] = box
        self.patch_set[box.id_box] = plt.Rectangle(
            (box.x, box.y),
            box.diam_x, box.diam_y,
            fill=False,
            edgecolor='black',
            linestyle='solid',
            lw=0.1
        )

    def get_patches(self):
        patch_list = []
        for x in self.patch_set.values():
            patch_list.append(x)
        return patch_list

    def get_box(self, id_box):
        return self.box_set[id_box], self.patch_set[id_box]

    def remove_box(self, id_box):
        if self.pending_box is not None:
            del self.box_set[self.pending_box]
            del self.patch_set[self.pending_box]
            self.pending_box = id_box
        else:
            self.pending_box = id_box


def ibex(q):
    proc = subprocess.Popen([
        "./optimizer-mop", "../benchs/MOP/binh.txt", "-f", "hc4", "-b",
        "largestfirst", "-s", "NDSdist", "--nb_ub_sols=3",
        "--linear-relax=compo", "--time=3600", "--eps_rel=0.01", "--w2=0.01"
        ], stdin=None, stdout=subprocess.PIPE
    )
    # box_set = {}
    state = "run"
    while state == "run":
        for line in proc.stdout:
            if 'add' in line.decode() or 'del' in line.decode():
                q.put(line.decode())
        proc.stdout.close()
        proc.wait()
        break


def plot():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_aspect(1)
    line, = plt.plot([], [], 'ro')
    line2, = plt.plot([], [], 'r-')
    plt.setp(line, color='b', linewidth=2.0, label='UB', marker='.')
    # line.set_data([], [])
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)
    return fig, line, ax, line2


def updateplot(num, q, l, data, ax, patch_list, l2):
    if pause:
        patch_list = [ax.add_patch(x) for x in data.get_patches()]
        return tuple(patch_list) + (l, ) + (l2, )
    else:
        try:
            result = q.get_nowait()

            if result != 'Q':
                # ax.clear()
                print(result)
                if 'add:' in result:
                    box = Box(result.split("add: ")[1])
                    data.append_box(box)
                if 'del:' in result:
                    id_box = eval(result.split("del: ")[1])['id']
                    data.remove_box(id_box)
                if 'add ub:' in result:
                    ub = UB(result.split("add ub: ")[1])
                    data.append_ub(ub)
                if 'del ub:' in result:
                    key = result.split("del ub: ")[1]
                    data.remove_ub(key)
                if 'add lb:' in result:
                    print(result)
                    lb = LB(result.split("add lb: ")[1])
                    data.append_lb(lb)
                l.set_data(data.get_ub_x(), data.get_ub_y())
                patch_list = [ax.add_patch(x) for x in data.get_patches()]
                l2.set_data(data.get_lb())
                return tuple(patch_list) + (l, ) + (l2, )
            else:
                print('done')
        except Exception as e:
            print(e)
            return l,


pause = False


def onClick(event, line):
    global pause
    if(event.button == 2):
        plt.xlim(min(line.get_xdata())-min(line.get_xdata())*.2,
                 max(line.get_xdata())*1.2)
        plt.ylim(min(line.get_ydata())-min(line.get_ydata())*.2,
                 max(line.get_ydata())*1.2)
    if(event.button == 3):
        pause ^= True


def main():
    q = multiprocessing.Queue()
    ibex_proc = multiprocessing.Process(None, ibex, args=(q,))
    ibex_proc.start()
    fig, line, ax, line2 = plot()
    patch_list = []
    fig.canvas.mpl_connect('button_press_event', lambda event: onClick(event, line))
    data = DataDict()
    line_ani = animation.FuncAnimation(
        fig, updateplot, q.qsize(), fargs=(q, line, data, ax, patch_list, line2),
        interval=1, blit=True, repeat=False
    )
    plt.show()
    print("Done")


if __name__ == '__main__':
    main()
