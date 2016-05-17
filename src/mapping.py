# Jialun Bao
# 11/13/15

import numpy as num
import matplotlib.pyplot as plt
import numpy.linalg as LA
import time as time

# constants
p = num.pi
A = 0.1
epsilon = 0.25
w = 0.6*p
Delta = 0.00001
delta = 0.05
dt = 0.05
T = 15
L = 400
H = 200
St = 0.05
R = 0
St_inverse = 1/St


start_time = time.time()

# wave function that defines the characteristics of 
# double gyre
def phi(x,y,t):
    temp = A*num.sin(p*f(x,t))*num.sin(p*y)
    return temp

def f(x,t):
    temp = epsilon*num.sin(w*t)*x**2+(1-2*epsilon*num.sin(w*t))*x
    return temp

# function that computes velocity of flow at each point
def velocity(r,t):
    x = r[:,0]
    y = r[:,1]
    vx = (phi(x,y+Delta,t)-phi(x,y-Delta,t))/(2*Delta)
    vy = (phi(x+Delta,y,t)-phi(x-Delta,y,t))/(2*Delta)
    return num.column_stack((-1*vx,vy))

def u_acc(r,t):
    return (velocity(r,t + Delta)- velocity(r,t- Delta))/(2*Delta)

def accel(r,v,t):
        a = (velocity(r,t) - v)*St_inverse + 1.5*R*u_acc(r,t)
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
        R[:,0] = num.clip(R[:,0],0.001,1.999)
        R[:,1] = num.clip(R[:,1],0.001,0.999)
        return num.hstack((R,V))
                

        

# y-coordinate of the grid
h = num.linspace(0.01,0.99,H)
# evolution time 0-10s
#tau = num.linspace(0,9.9,100)

# getting 100 mapping data files
for n in range(0,1):
    #output = open('mapping%d.txt'%n,'ab')
    output = open('R0_St005.txt','ab')
    for i in h:
        
        # initialize grid horizontally
        x = num.linspace(0.01,1.99,L)
        y = num.linspace(i,i,L)
        r = num.column_stack((x,y))
        #print(r)
        #v = num.zeros((L,2))
        state = num.hstack((r,velocity(r,0)))
        #print(state)
        
        # perform RK4 to get position of particle 20s later
        for t in num.arange(0,T,dt):
            state = rk4(state,t)
                         
        # append data to the file       
        num.savetxt(output,state[:,0:2])
    output.close()

print(time.time()-start_time)



