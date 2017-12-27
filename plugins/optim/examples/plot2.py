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
import traceback


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __lt__(self, other):
        if self.x < other.x:
            return True
        elif self.x == other.x:
            return self.y < other.y
        return False


class Box:
    def __init__(self, box_string):
        box_string = box_string.replace('inf', "math.inf")
        box_string = box_string.replace('nan', "0")
        box_var = eval(box_string)
        self.createBox(box_var['id'], box_var['pts'][0], box_var['pts'][1],
                       box_var['diam_x'], box_var['diam_y'], box_var['pA'],
                       box_var['pB'])

    def __repr__(self):
        return str((self.x, self.y))

    def __lt__(self, other):
        if self.x < other.x:
            return True
        elif self.x == other.x:
            return self.y < other.y
        return False

    def createBox(self, id_box, x, y,  diam_x, diam_y, pA, pB):
        self.id_box = id_box
        self.x = x
        self.y = y
        self.diam_x = diam_x
        self.diam_y = diam_y
        self.pA = pA
        self.pB = pB


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

    def get_ub(self):
        ub_list_x = []
        ub_list_y = []
        ub_list = []
        for x in self.UB_set.values():
            ub_list.append(Point(x.x, x.y))
        ub_list = sorted(ub_list)
        i = 0
        if ub_list:
            ub_list_x.append(ub_list[0].x)
            ub_list_y.append(ub_list[0].y)
            for point in ub_list[1:]:
                ub_list_x.append(point.x)
                ub_list_y.append(ub_list[i].y)
                ub_list_x.append(point.x)
                ub_list_y.append(point.y)
                i = i + 1
            box = sorted(list(self.box_set.values()))[-1]
            ub_list_x.append(box.x + box.diam_x)
            ub_list_y.append(ub_list_y[-1])
        return ub_list_x, ub_list_y

    def get_lb(self):
        lb_list_x = []
        lb_list_y = []
        box_list = []
        for x in self.LB_set.values():
            box_list.append(Point(x.x, x.y))
        for x in self.box_set.values():
            box_list.append(Point(x.x, x.y))
        box_list = sorted(box_list)
        i = 0
        for box in box_list[1:]:
            if box.y > box_list[i].y:
                box_list.remove(box)
            else:
                i = i + 1
        i = 0
        if box_list:
            lb_list_x.append(box_list[0].x)
            lb_list_y.append(box_list[0].y)
            for box in box_list[1:]:
                if box.x != box_list[i].x:
                    lb_list_x.append(box.x)
                    lb_list_y.append(box_list[i].y)
                lb_list_x.append(box.x)
                lb_list_y.append(box.y)
                i = i + 1
            box = sorted(list(self.box_set.values()))[-1]
            lb_list_x.append(box.x + box.diam_x)
            lb_list_y.append(lb_list_y[-1])
        return lb_list_x, lb_list_y

    def get_cy_lines(self):
        line_list = []
        for box in self.box_set.values():
            if(box.pA[0] < box.pB[0]) and (box.pA[1] > box.pB[1]):
                line = plt.Line2D((box.pA[0], box.pB[0]),
                                  (box.pA[1], box.pB[1]),
                                  lw=0.5, markeredgecolor='black')
                line_list.append(line)
        return line_list

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
            lw=0.1,
            label="Box"
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
        "./optimizer-mop", "../benchs/MOP/osy.txt", "-f", "acidhc4", "-b",
        "largestfirst", "-s", "NDSdist", "--nb_ub_sols=50",
        "--linear-relax=compo", "--time=3600", "--eps=1", "--w2=0.01",
        "--cy-contract-full"
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
    line, = plt.plot([], [], 'b-', label='Upperbound')
    line2, = plt.plot([], [], 'r-', label='Lowerbound')
    plt.xlabel('Funcion Objetivo 1')
    plt.ylabel('Funcion Objetivo 2')
    plt.title('HECHO POR MATIAS CAMPUSANO')
    # plt.setp(line, color='b', linewidth=2.0, label='UB', marker='-')
    # line.set_data([], [])
    plt.xlim(-300, 300)
    plt.ylim(-300, 300)
    return fig, line, ax, line2


def updateplot(num, q, l, data, ax, patch_list, l2, line_list):
    if pause:
        patch_list = [ax.add_patch(x) for x in data.get_patches()]
        return tuple(patch_list) + (l, l2,) + tuple(line_list)
    else:
        try:
            result = q.get_nowait()
            if result != 'Q':
                # ax.clear()
                # print(result)
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
                    lb = LB(result.split("add lb: ")[1])
                    data.append_lb(lb)
                l.set_data(data.get_ub())
                patch_list = [ax.add_patch(x) for x in data.get_patches()]
                line_list = [ax.add_line(x) for x in data.get_cy_lines()]
                l2.set_data(data.get_lb())
                # plt.legend(handles=[l, l2, patch_list[0]])
                return tuple(patch_list) + (l, ) + (l2, ) + tuple(line_list)
            else:
                print('done')
        except Exception as e:
            print("Esto es un error:")
            print(e)
            return tuple(patch_list) + (l, ) + (l2, ) + tuple(line_list)


pause = False


def onClick(event, line, line2):
    global pause
    if(event.button == 2):
        if min(line.get_xdata()[:-1]) < min(line2.get_xdata()):
            x_min = min(line.get_xdata())
        else:
            x_min = min(line2.get_xdata())
        if max(line.get_xdata()[:-1]) > max(line2.get_xdata()):
            x_max = max(line.get_xdata())
        else:
            x_max = max(line2.get_xdata())
        if x_min < 0:
            x_min = x_min*1.5
        else:
            x_min = x_min*.8
        if x_max < 0:
            x_max = x_max*0.8
        else:
            x_max = x_max*1.5
        if min(line2.get_ydata()) < 0:
            y_min = min(line2.get_ydata())*1.5
        else:
            y_min = min(line2.get_ydata())*.8
        if max(line.get_ydata()) < 0:
            y_max = max(line.get_ydata())*0.8
        else:
            y_max = max(line.get_ydata())*1.5
        plt.xlim(x_min,
                 x_max)
        plt.ylim(y_min,
                 y_max)
    if(event.button == 3):
        pause ^= True


def main():
    try:
        q = multiprocessing.Queue()
        ibex_proc = multiprocessing.Process(None, ibex, args=(q,))
        ibex_proc.start()
        fig, line, ax, line2 = plot()
        patch_list = []
        line_list = []
        fig.canvas.mpl_connect('button_press_event',
                               lambda event: onClick(event, line, line2))
        data = DataDict()
        line_ani = animation.FuncAnimation(
            fig, updateplot, q.qsize(),
            fargs=(q, line, data, ax, patch_list, line2, line_list),
            interval=0, blit=True, repeat=False
        )
        plt.show()
        print("Done")
        ibex_proc.terminate()
    except Exception:
        ibex_proc.terminate()


if __name__ == '__main__':
    main()
