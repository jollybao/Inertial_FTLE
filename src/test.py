import pylab as plt
import numpy as num
import matplotlib.animation as animation


fig, ax = plt.subplots(1,1,figsize=(10,5))
#plt.rcParams['animation.ffmpeg_path'] = 'C:/ffmpeg/bin/ffmpeg'
#mywriter = animation.FFMpegWriter()

# constants
p = num.pi
A = 0.1
epsilon = 0.25
w = 0.6*p
delta = 0.0001
dt = 0.1
partition = 20
St = 1
R = 1
St_inverse = 1/St
#N = int(input("Enter number of particles: "))
col = ['r','y','b','g','k','c','m','r','y','b',
       'g','k','c','m','r','y','b','g','k','c','m',
       'r','y','b','g','k','c','m']

# wave function that defines the characteristics of 
# double gyre
def phi(x,y,t):
    temp = A*num.sin(p*f(x,t))*num.sin(p*y)
    return temp

def f(x,t):
    temp = epsilon*num.sin(w*t)*x**2+(1-2*epsilon*num.sin(w*t))*x
    return temp

def velocity(x,y,t):
    vx = (phi(x,y+delta,t)-phi(x,y-delta,t))/(2*delta)
    vy = (phi(x-delta,y,t)-phi(x+delta,y,t))/(2*delta)
    return -1*vx,-1*vy

def vel(r,t):
    x = r[:,0]
    y = r[:,1]
    vx = (phi(x,y+delta,t)-phi(x,y-delta,t))/(2*delta)
    vy = (phi(x-delta,y,t)-phi(x+delta,y,t))/(2*delta)
    return num.column_stack((-1*vx,-1*vy))


def accel(x,y,t):
    vx1 = (phi(x,y+delta,t+delta)-phi(x,y-delta,t+delta))/(2*delta)
    vx2 = (phi(x,y+delta,t-delta)-phi(x,y-delta,t-delta))/(2*delta)
    vy1 = (phi(x-delta,y,t+delta)-phi(x+delta,y,t+delta))/(2*delta)
    vy2 = (phi(x-delta,y,t-delta)-phi(x+delta,y,t-delta))/(2*delta)
    ax = (vx1-vx2)/(2*delta)
    ay = (vy1-vy2)/(2*delta)
    ux,uy = velocity(x,y,t)
    ux = ux*St_inverse
    uy = uy*St_inverse
    return ux-1.5*ax,uy-1.5*ay    


    
# make a 2D mesh grid of size 40*20
X,Y = plt.meshgrid(num.arange(0,2,1/partition),num.arange(0,1,1/partition))
Vx,Vy = velocity(X,Y,0.1)

Ax,Ay = accel(X,Y,0)

# vector arrows
Q = ax.quiver(X,Y,Ax,Ay,scale=10)
plt.show()

