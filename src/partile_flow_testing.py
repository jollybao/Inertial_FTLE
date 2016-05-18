# Jialun Bao
# 11/13/15

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
delta = 0.00001
dt = 0.1
partition = 20
St = 1
R = 1
St_inverse = 1/St
N = int(input("Enter number of particles: "))
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

def u_acc(r,t):
    return (vel(r,t + delta)- vel(r,t- delta))/(2*delta)

def accel(r,v,t):
        a = 1.5*R*u_acc(r,t)
        #(vel(r,t) - v)*St_inverse +
        return a

def rk4(state,t):
        r = state[:,0:2]
        v = state[:,2:4]
        a = accel(r,v,t)
            
        r1 = r + 0.5*v*dt
        v1 = v + 0.5*a*dt
        a1 = accel(r1,v1,t+dt/2)
        
        r2 = r + 0.5*v1*dt
        v2 = v + 0.5*a1*dt
        a2 = accel(r2,v2,t+dt/2)
        
        r3 = r + v2*dt
        v3 = v + a2*dt
        a3 = accel(r3,v3,t+dt)
        
        R = r + (dt/6)*(v + 2*v1 + 2*v2 + v3)
        V = v + (dt/6)*(a + 2*a1 + 2*a2 + a3)
        #R[:,0] = num.clip(R[:,0],0,2)
        #R[:,1] = num.clip(R[:,1],0,1)
        return num.hstack((R,V))
    
def NonInertial_rk4(state,t):
        r = state[:,0:2]
        k1 = dt*vel(r,t)
        k2 = dt*vel(r+0.5*k1,t+0.5*dt)
        k3 = dt*vel(r+0.5*k2,t+0.5*dt)
        k4 = dt*vel(r+k3,t+dt)
        r += (k1+2*k2+2*k3+k4)/6
        state[:,0:2] = r
        return state

    
# make a 2D mesh grid of size 40*20
X,Y = plt.meshgrid(num.arange(0,2,1/partition),num.arange(0,1,1/partition))
Vx,Vy = velocity(X,Y,0.1)

# vector arrows
Q = ax.quiver(X,Y,Vx,Vy,scale=10)

# initialize array of particles
C = num.empty([N],plt.Circle)
for i in range(0,N):
    C[i] = plt.Circle((-1,-1),radius = 0.03, fc = col[i])
    
state = num.empty([N,4],float)
for i in range(0,N):
    #print("Enter x and y coordinates of the circle ",i+1)
    state[i,0] = (i*0.5 + 0.1)%2
    state[i,1] = (i//4)*0.3 + 0.1
    state[i,2],state[i,3] = velocity(state[i,0],state[i,1],0)
    
    C[i].center = (state[i][0],state[i][1])
    ax.add_patch(C[i])
    

# animation for particle moving along the vector field
def animate(num,Q,X,Y,C,state,N):
    t = num/1
    dt = 1/10
    Vx,Vy = velocity(X,Y,t)
    Q.set_UVC(Vx,Vy)  
    
    # update particles' positions
    for j in range(0,10):
        state = rk4(state,t)
            
    for i in range(0,N):
        C[i].center = (state[i][0],state[i][1])
    return Q,C

ani = animation.FuncAnimation(fig, animate,
         fargs=(Q,X,Y,C,state,N),
    interval=100,blit=False)

#ani.save('VF_demo.mp4',writer = mywriter)

plt.show()

