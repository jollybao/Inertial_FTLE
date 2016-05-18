import numpy as num
import matplotlib.pyplot as plt
import numpy.linalg as LA
import time as time

# constants
L = 400
H = 200
delta_1 = 1.98/(L-1)
delta_2 = 0.98/(H-1)

start_time = time.time()

# spatial Jacobian that is used to compute FTLE
def Jacobian(X,Y):
    J = num.empty([2,2],float)
    FTLE = num.empty([H-2,L-2],float)
    
    for i in range(0,H-2):
        for j in range(0,L-2):
            J[0][0] = (X[(1+i)*L+2+j]-X[(1+i)*L+j])/(2*delta_1)
            J[0][1] = (X[(2+i)*L+1+j]-X[i*L+1+j])/(2*delta_1)
            J[1][0] = (Y[(1+i)*L+2+j]-Y[(1+i)*L+j])/(2*delta_2)
            J[1][1] = (Y[(2+i)*L+1+j]-Y[i*L+1+j])/(2*delta_2)
			
			# Green-Cauchy tensor
            D = num.dot(num.transpose(J),J)
			# its largest eigenvalue
            lamda = LA.eigvals(D)
            FTLE[i][j] = max(lamda)
    return FTLE


start_time = time.time()
'''
#Input = open('heavy_particle.txt','r')
Input = open('test10.txt','r')
X,Y = num.loadtxt(Input,unpack = True)
Input.close()

plt.figure(1)
plt.scatter(X,Y)
'''

input1 = open('rk4.txt','r')
input2 = open('rk4_accel.txt','r')
x1,y1,vx1,vy1 = num.loadtxt(input1,unpack = True)
x2,y2,vx2,vy2 = num.loadtxt(input2,unpack = True)
input1.close()
input2.close()

x1 = x1[::2]
y1 = y1[::2]
vx1 = vx1[::2]
vy1 = vy1[::2]
x2 = x2[::2]
y2 = y2[::2]
vx2 = vx2[::2]
vy2 = vy2[::2]
plt.figure(1)
plt.scatter(x1,y1)

plt.figure(2)
plt.scatter(x2,y2)


plt.figure(3)
plt.plot(vx1,'k')
plt.plot(vx2,'b')

plt.figure(4)
plt.plot(vy1,'k')
plt.plot(vy2,'b')

#f, (ax1, ax2) = plt.subplots(2, 1, sharey=True)
#ax1.scatter(x1, y1)

#ax2.scatter(x2, y2)
'''
plt.figure(2)
FTLE = Jacobian(X,Y)
FTLE = num.log(FTLE)
plt.imshow(FTLE)
plt.colorbar()
plt.gca().invert_yaxis()
'''
plt.show()

print(time.time()-start_time)
plt.show()
