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
w = 0.2*p
Delta = 0.00001
dt = 0.1
delta = dt
T = 15
L = 400
H = 200
St = 1
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
def velocity(x,y,t):
    f = epsilon*num.sin(w*t)*x*x + (1-2*epsilon*num.sin(w*t))*x
    vx = A*p*num.sin(p*f)*num.cos(p*y)
    vy = -1*A*p*num.cos(p*f)*num.sin(p*y)*(2*epsilon*num.sin(w*t)*x + (1-2*epsilon*num.sin(w*t)))
    
    return num.column_stack((vx,vy))

def u_acc(x,y,t):
    u_accel = (velocity(x,y,t + Delta) - velocity(x,y,t - Delta))/(2*Delta)
    v = velocity(x,y,t)
    
    u_ax = (velocity(x+Delta,y,t) - velocity(x-Delta,y,t))/(2*Delta)
    u_ay = (velocity(x,y+Delta,t) - velocity(x,y-Delta,t))/(2*Delta)
    u_accel = u_accel + u_ax*num.column_stack((v[:,0],v[:,0])) + u_ay*num.column_stack((v[:,1],v[:,1]))
    return u_accel


def accel(x,y,v,t):
    u = velocity(x,y,t)
    a = St_inverse*(u - v)
    au = 1.5*R*u_acc(x,y,t)
    a = a + au
    return a

def update(state,t):
    r = state[:,2:4]
    v = accel(state[:,0],state[:,1],r,t)
    return num.hstack((r,v))

    
def rk4(state,t):
    r = state
    k1 = dt*update(r,t)
    k2 = dt*update(r+0.5*k1,t+0.5*dt)
    k3 = dt*update(r+0.5*k2,t+0.5*dt)
    k4 = dt*update(r+k3,t+dt)
    r += (k1+2*k2+2*k3+k4)/6
    state[:,0] = num.clip(r[:,0],0.00001,1.99999)
    state[:,1] = num.clip(r[:,1],0.00001,0.99999)
    return state

'''
r = num.array([[0.9,0.1],[1,0.4]],num.float64)
state = num.hstack((r,velocity(r[:,0],r[:,1],0)))
pts = num.zeros((150,2))
i = 0

for t in num.arange(0,T,dt):
    pts[i,:] = state[0,0:2]
    i = i + 1
    state = test(state,t)


plt.figure(1)
plt.plot(pts[:,0],pts[:,1])

plt.show()
'''

output = open('R1_St1.txt','ab')
# y-coordinate of the grid
h = num.linspace(0.01,0.99,H,num.float64)
#pts = num.zeros((H*L,2))
#j = 0
for i in h:
    
    # initialize grid horizontally
    x = num.linspace(0.01,1.99,L,num.float64)
    y = i*num.ones(L,num.float64)
    r = num.column_stack((x,y))
    state = num.hstack((r,velocity(x,y,0)))
    #print(state)
    
    # perform RK4 to get position of particle 20s later
    for t in num.arange(0,T,dt):
  
        state = rk4(state,t)
        #state = NonInertial_rk4(state,t)
    #pts[j:L+j,:] = state[:,0:2]
    #j = j + L
    num.savetxt(output,state[:,0:2])
       
#plt.scatter(pts[:,0],pts[:,1])
#plt.show()    


output.close()
print(time.time()-start_time)



