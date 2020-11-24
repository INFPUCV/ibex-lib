import paramiko
import getpass
import numpy as np

print("Cargando ibex...")

home = "/home/iaraya/ibex/ibex-lib"

ssh = None
def connect(host):
    global ssh
    port = 22
    username = input("login:")
    password = getpass.getpass("pass:")

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host, port, username, password)


def search_efficient(instance, box_y, mode="efficient", prec=1e-5, box_x=None):
    boxy_str = [str(y) for y in box_y]
    if box_x:
        boxx_str = [str(x) for x in box_x]
    
    if box_x==None:
        stdin, stdout, stderr = ssh.exec_command(home+"/__build__/plugins/optim-mop/search_efficient "+ home + \
                "/plugins/optim-mop/benchs/" + instance + " --eps=" + str(prec) + \
                " --box=\"" + " ".join(boxy_str) +"\" --se_mode=" + mode)
    else:
         stdin, stdout, stderr = ssh.exec_command(home+"/__build__/plugins/optim-mop/search_efficient "+ home + \
                "/plugins/optim-mop/benchs/" + instance + " --eps=" + str(prec) + \
                " --box=\"" + " ".join(boxy_str) + "\" --box_x=\"" + " ".join(boxx_str) + "\" --se_mode=" + mode)
    lines = stdout.readlines()
    values_str = lines[0].split()
    values = [float(v) for v in values_str]
    y = np.array(values[0:2])
    x = np.array(values[2:len(values)])
    time = float(lines[1].strip())
    
    #TODO: verificar que "y" sea una solución eficiente 
    
    return x, y, time

#llama al solver para encontrar frontera de pareto
def solve(instance, prec=1e-1):
    cmd = home+"/__build__/plugins/optim-mop/ibexmop "+ home + \
                "/plugins/optim-mop/benchs/" + instance + " --eps=" + str(prec) + \
                " --eps-contract --cy-contract-full --ub=ub1 --verbose"

    #print(cmd) 
    stdin, stdout, stderr = ssh.exec_command(cmd)
    
    flag = False; y1 = None; y2 = None
    i = 0
    
    for line in stdout.readlines():
        if "number of solutions:" in line:
            aux, n = line.strip().split(':')
            sols = int(n)
            y1 = np.zeros(sols); y2 = np.zeros(sols)
            
        if line==" solutions:\n":
            flag=True; continue
        if flag:
            a,b = line.strip().split()
            y1[i] = float(a); y2[i] = float(b)
            i += 1
    
    return y1, y2


def get_instances():
    cmd = "ls "+home+"/plugins/optim-mop/benchs/*.txt  | xargs -n 1 basename"
    stdin, stdout, stderr = ssh.exec_command(cmd)
    for line in stdout.readlines():
        print(line.strip())
        

## MOP-SERVER ###
def init_mopserver(instance, port=8000):
    cmd = "killall ibexmop-server; "+home+"/__build__/plugins/optim-mop/ibexmop-server "+home+"/plugins/optim-mop/benchs/"+instance+ \
                     " --cy-contract-full --port="+str(port)+" --server_mode --ub=ub1"
    print(cmd)
    transport = ssh.get_transport()
    channel = transport.open_session()
    channel.exec_command(cmd)
    return channel

def close_mopserver(port=8000):
    stdin, stdout, stderr = ssh.exec_command("echo fns | netcat localhost "+str(port))
    print(stdout.readlines())

def print_lines(stdout):
    for line in stdout.readlines():
        print(line)
    
def run(iters, prec=1e-2, port=8000):
    print("echo run "+str(iters)+" "+str(prec)+" | netcat localhost "+str(port))
    stdin, stdout, stderr = ssh.exec_command("echo run "+str(iters)+" "+str(prec)+" | netcat localhost "+str(port))
    print_lines(stdout)
    
    
def update_refpoint(y1, y2, prec=1e-8, port=8000):
    stdin, stdout, stderr = ssh.exec_command("echo zoo "+str(y1)+ " "+str(y2)+" "+ str(prec) +" | netcat localhost "+str(port))
    print_lines(stdout)
    
def get_envelope(prec=1e-8, port=8000, l1=-1e8, l2=1e8, u1=-1e8, u2=1e8):
    u1 = []; u2 = []
    stdin, stdout, stderr = ssh.exec_command("echo gup 1e-8 -1e8 1e8 -1e8 1e8 | netcat localhost "+str(port))   
    for point in stdout.readlines()[0].split('/'):
        if (len(point.split(',')) == 2):
            uu1 ,uu2 = point.split(',')
            u1.append(float(uu1))
            u2.append(float(uu2))
            
    l1 = []; l2 = []
    stdin, stdout, stderr = ssh.exec_command("echo glw 1e-8 -1e8 1e8 -1e8 1e8 | netcat localhost "+str(port))
    for point in stdout.readlines()[0].split('/'):
        if (len(point.split(',')) == 2):
            ll1 ,ll2 = point.split(',')
            l1.append(float(ll1))
            l2.append(float(ll2))
    
    return np.array(l1),np.array(l2), np.array(u1),np.array(u2) 

import matplotlib.patches as patches

def get_boxes(port=8000):
    stdin, stdout, stderr = ssh.exec_command("echo get_boxes | netcat localhost "+str(port))
    boxes = []
    for line in stdout.readlines() :
        r = [float(x) for x in line.split(' ')]
        r = np.array(r)
        boxes.append(r)
        
    return np.array(boxes)

def add_boxes(boxes, plt):
    # Create figure and axes
    fig,ax = plt.subplots(1)
    for r in boxes:
        rect = patches.Rectangle( ( r[0], r[2]) ,(r[1]-r[0]),(r[3]-r[2]), linewidth=1,edgecolor='r', facecolor=None, alpha=0.1) # ,facecolor='r', alpha=0.1)
        # Add the patch to the Axes
        ax.add_patch(rect)