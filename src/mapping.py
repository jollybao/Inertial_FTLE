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
Delta = 0.000001
delta = 0.005
dt = 0.005
T = 15
L = 40
H = 20
St = 10
R = 1
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
        a = R*u_acc(r,t)
        #(velocity(r,t) - v)*St_inverse
        return a

def update(state,t):
    r = state[:,2:4]
    v = accel(state[:,0:2],r,t)
    return num.hstack((r,v))

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
        R[:,0] = num.clip(R[:,0],0.01,1.99)
        R[:,1] = num.clip(R[:,1],0.01,0.99)
        return num.hstack((R,V))
                
def NonInertial_rk4(state,t):
        r = state[:,0:2]
        k1 = 1.5*dt*velocity(r,t)
        k2 = 1.5*dt*velocity(r+0.5*k1,t+0.5*dt)
        k3 = 1.5*dt*velocity(r+0.5*k2,t+0.5*dt)
        k4 = 1.5*dt*velocity(r+k3,t+dt)
        r += (k1+2*k2+2*k3+k4)/6
        state[:,0:2] = r
        state[:,2:4] = velocity(r,t+dt)
        return state
    
def test(state,t):
    k1 = dt*update(state,t)
    k2 = dt*update(state+0.5*k1,t+0.5*dt)
    k3 = dt*update(state+0.5*k2,t+0.5*dt)
    k4 = dt*update(state+k3,t+dt)
    state += (k1+2*k2+2*k3+k4)/6
    return state

def test2(state,t):
    r = state[:,0:2]
    v = state[:,2:4]

    k1 = dt*accel(r,v,t)
    l1 = dt*v

    k2 = dt*accel(r+l1/2,v+k1/2,t+dt/2)
    l2 = dt*(v+k1/2)

    k3 = dt*accel(r+l2/2,v+k2/2,t+dt/2)
    l3 = dt*(v+k2/2)
    
    k4 = dt*accel(r+l3,v+k3,t+dt)
    l4 = dt*(v+k3)
    
    

    R = r + (l1 + 2*l2 + 2*l3 + l4)/6
    V = v + (k1 + 2*k2 + 2*k3 + k4)/6
    R[:,0] = num.clip(R[:,0],0.001,1.999)
    R[:,1] = num.clip(R[:,1],0.001,0.999)
    return num.hstack((R,V))    
    
# y-coordinate of the grid
h = num.linspace(0.01,0.99,H,num.float64)
# evolution time 0-10s
#tau = num.linspace(0,9.9,100)

'''
# getting 100 mapping data files
for n in range(0,1):
    #output = open('mapping%d.txt'%n,'ab')
    output = open('test10.txt','ab')
    for i in h:
        
        # initialize grid horizontally
        x = num.linspace(0.01,1.99,L,num.float64)
        y = num.linspace(i,i,L,num.float64)
        r = num.column_stack((x,y))
        #print(r)
        #v = num.zeros((L,2))
        state = num.hstack((r,velocity(r,0)))
        #print(accel(state[:,0:2],state[:,2:4],0))
        #print(state)
        
        # perform RK4 to get position of particle 20s later
        for t in num.arange(0,T,dt):
            state = test2(state,t)
                         
        # append data to the file       
        num.savetxt(output,state[:,0:2])
    output.close()
'''

output1 = open('rk4_accel.txt','ab')
#output2 = open('rk4.txt','ab')
r = num.array([[1,0.5],[1,0.4]],num.float64)
state1 = num.hstack((r,velocity(r,0)))
state2 = state1
for t in num.arange(0,T,dt):
            #print(state1[0,:])
            num.savetxt(output1,state1)
            state1 = test2(state1,t)
            
            #print(state2[0,:])
            #num.savetxt(output2,state2)
            #state2 = NonInertial_rk4(state2,t)
            


num.savetxt(output1,state1)
#num.savetxt(output2,state2)
output1.close()
#output2.close()
print(time.time()-start_time)



