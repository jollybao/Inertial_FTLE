import numpy as num
import matplotlib.pyplot as plt
import numpy.linalg as LA
import time as time

# constants
L = 400
H = 200
# horizontal grid difference
delta_1 = 1.98/(L-1)
# verticle grid difference
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


Input = open('step001.txt','r')
X,Y = num.loadtxt(Input,unpack = True)
Input.close()

# scatter plot
plt.figure(1)
plt.scatter(X,Y)


# FTLE plot
plt.figure(2)
FTLE = Jacobian(X,Y)
FTLE = num.log(FTLE)
plt.imshow(FTLE)
plt.colorbar()
#plt.gca().invert_yaxis()


print(time.time()-start_time)
plt.show()
