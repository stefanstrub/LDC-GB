import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
from mpl_toolkits import mplot3d
from copy import deepcopy

def Gradient_Descent(Grad,x,y, gamma = 0.00125, epsilon=0.0001, nMax = 10000 ):
    #Initialization
    i = 0
    iter_x, iter_y, iter_count = np.empty(0),np.empty(0), np.empty(0)
    error = 10
    X = np.array([x,y])
    
    #Looping as long as error is greater than epsilon
    while np.linalg.norm(error) > epsilon and i < nMax:
        i +=1
        iter_x = np.append(iter_x,x)
        iter_y = np.append(iter_y,y)
        iter_count = np.append(iter_count ,i)   
        #print(X) 
        
        X_prev = X
        X = X - gamma * Grad(x,y)
        error = X - X_prev
        x,y = X[0], X[1]
          
    print(X)
    return X, iter_x,iter_y, iter_count

def Coordinate_Descent(Grad,x,y, gamma = 0.00125, epsilon=0.0001, nMax = 10000 ):
    #Initialization
    i = 0
    iter_x, iter_y, iter_count = np.empty(0),np.empty(0), np.empty(0)
    error = 10
    X = np.array([x,y])
    
    #Looping as long as error is greater than epsilon
    while np.linalg.norm(error) > epsilon and i < nMax:
        i +=1
        iter_x = np.append(iter_x,x)
        iter_y = np.append(iter_y,y)
        iter_count = np.append(iter_count ,i)   
        #print(X) 
        
        X_prev = deepcopy(X)
        print(X) 
        x = X[0] - gamma * Grad(x,y)[0]
        y = X[1] - gamma * Grad(x,y)[1]
        X = np.array([x,y])
        error = X - X_prev
        # x,y = X[0], X[1]
          
    print('CD',X, error)
    return X, iter_x,iter_y, iter_count

def Coordinate_Monte_Carlo(function,x,y, gamma = 0.00125, epsilon=0.0001, nMax = 10000 ):
    #Initialization
    i = 0
    iter_x, iter_y, iter_count = np.empty(0),np.empty(0), np.empty(0)
    error = 10
    X = np.array([x,y])
    
    #Looping as long as error is greater than epsilon
    while np.linalg.norm(error) > epsilon and i < nMax:
        i +=1
        iter_x = np.append(iter_x,x)
        iter_y = np.append(iter_y,y)
        iter_count = np.append(iter_count ,i)   
        #print(X) 
        
        X_prev = deepcopy(X)
        size = 10
        suggesstions = np.random.uniform(low=-5,high=5,size=size)
        suggesstions[0] = X[0]
        suggesstions = np.transpose([suggesstions,np.ones(size)*X[1]])
        results = np.zeros(size)
        for n in range(len(suggesstions)):
            results[n] = function(*suggesstions[n])
        X = suggesstions[np.argmin(results)]

        suggesstions = np.random.uniform(low=-5,high=5,size=size)
        suggesstions[1] = X[1]
        suggesstions = np.transpose([np.ones(size)*X[0],suggesstions])
        results = np.zeros(size)
        for n in range(len(suggesstions)):
            results[n] = function(*suggesstions[n])
        X = suggesstions[np.argmin(results)]
        if np.linalg.norm(X - X_prev) != 0:
            error = X - X_prev
        x,y = X[0], X[1]
          
    print('CMC',X, error)
    return X, iter_x,iter_y, iter_count


def Newton_Raphson_Optimize(Grad, Hess, x,y, epsilon=0.000001, nMax = 200):
    #Initialization
    i = 0
    iter_x, iter_y, iter_count = np.empty(0),np.empty(0), np.empty(0)
    error = 10
    X = np.array([x,y])
    
    #Looping as long as error is greater than epsilon
    while np.linalg.norm(error) > epsilon and i < nMax:
        i +=1
        iter_x = np.append(iter_x,x)
        iter_y = np.append(iter_y,y)
        iter_count = np.append(iter_count ,i)   
        print(X) 
        
        X_prev = X
        X = X - np.linalg.inv(Hess(x,y)) @ Grad(x,y)
        error = X - X_prev
        x,y = X[0], X[1]
          
    return X, iter_x,iter_y, iter_count

def f_2(x,y):
    return .01*x**2 - .1*y**2

def Himmer(x,y):
    return (x**2 + y - 11)**2 + ( x + y**2 - 7 )**2

def Grad_Himmer(x,y):
    return np.array([2 * (-7 + x + y**2 + 2 * x * (-11 + x**2 + y)), 2 * (-11 + x**2 + y + 2 * y * (-7 + x + y**2))])

def Hessian_Himmer(x,y):
    h11 = 4 * (x**2 + y - 11) + 8 * x**2 + 2
    h12 = 4 * x + 4 * y
    h21 = 4 * x + 4 * y 
    h22 = 4 * (x + y**2 - 7) + 8 * y**2 + 2
    
    return np.array([[h11,h12],[h21,h22]]) 

start = [-2,-1]
root_gd,iter_x_gd,iter_y_gd, iter_count_gd = Gradient_Descent(Grad_Himmer,start[0],start[1],gamma = 0.001, epsilon=0.01, nMax = 1000)
print('Found min gd:', Himmer(iter_x_gd[-1],iter_y_gd[-1]), iter_count_gd[-1])
root_cd,iter_x_cd,iter_y_cd, iter_count_cd = Coordinate_Descent(Grad_Himmer,start[0],start[1],gamma = 0.001, epsilon=0.01, nMax = 1000)
print('Found min cd:', Himmer(iter_x_cd[-1],iter_y_cd[-1]), iter_count_cd[-1])
root_cmc,iter_x_cmc,iter_y_cmc, iter_count_cmc = Coordinate_Monte_Carlo(Himmer,start[0],start[1],gamma = 0.001, epsilon=0.01, nMax = 100)
print('Found min cmc:', Himmer(iter_x_cmc[-1],iter_y_cmc[-1]), iter_count_cmc[-1])
#(Grad,x,y, gamma = 0.00125, epsilon=0.0001, nMax = 10000 )

root_nr,iter_x_nr,iter_y_nr, iter_count_nr = Newton_Raphson_Optimize(Grad_Himmer,Hessian_Himmer,start[0],start[1], nMax = 50)
print('Found min nr:', Himmer(iter_x_nr[-1],iter_y_nr[-1]), iter_count_nr[-1])

x = np.linspace(-5,5,100)
y = np.linspace(-5,5,100)
X, Y = np.meshgrid(x, y)
Z = Himmer(X, Y)

#Angles needed for quiver plot
anglesx = iter_x_gd[1:] - iter_x_gd[:-1]
anglesy = iter_y_gd[1:] - iter_y_gd[:-1]
anglesx_cd = iter_x_cd[1:] - iter_x_cd[:-1]
anglesy_cd = iter_y_cd[1:] - iter_y_cd[:-1]
anglesx_cmc = iter_x_cmc[1:] - iter_x_cmc[:-1]
anglesy_cmc = iter_y_cmc[1:] - iter_y_cmc[:-1]
anglesx_nr = iter_x_nr[1:] - iter_x_nr[:-1]
anglesy_nr = iter_y_nr[1:] - iter_y_nr[:-1]

fig = plt.figure(figsize = (16,8))

#Surface plot
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_surface(X,Y,Z,rstride = 5, cstride = 5, cmap = 'jet', alpha = .4, edgecolor = 'none' )
ax.plot(iter_x_gd,iter_y_gd, f_2(iter_x_gd,iter_y_gd),color = 'orange', marker = '*', alpha = .4, label = 'Gradient descent')
ax.plot(iter_x_nr,iter_y_nr, f_2(iter_x_nr,iter_y_nr),color = 'darkblue', marker = 'o', alpha = .4, label = 'Newton')
ax.legend()

#Rotate the initialization to help viewing the graph
ax.view_init(45, 60)
ax.set_xlabel('x')
ax.set_ylabel('y')

#Contour plot
ax = fig.add_subplot(1, 2, 2)
ax.contour(X,Y,Z, 60, cmap = 'jet')

#Plotting the iterations and intermediate values
# ax.scatter(iter_x_gd,iter_y_gd,color = 'orange', marker = '*', label = 'Gradient descent')
# ax.quiver(iter_x_gd[:-1], iter_y_gd[:-1], anglesx, anglesy, scale_units = 'xy', angles = 'xy', scale = 1, color = 'orange', alpha = .3)
ax.scatter(iter_x_cd,iter_y_cd,color = 'green', marker = '*', label = 'Coordinate descent')
ax.quiver(iter_x_cd[:-1], iter_y_cd[:-1], anglesx_cd, anglesy_cd, scale_units = 'xy', angles = 'xy', scale = 1, color = 'orange', alpha = .3)
ax.scatter(iter_x_cmc,iter_y_cmc,color = 'red', marker = '*', label = 'Coordinate Monte Carlo')
ax.quiver(iter_x_cmc[:-1], iter_y_cmc[:-1], anglesx_cmc, anglesy_cmc, scale_units = 'xy', angles = 'xy', scale = 1, color = 'orange', alpha = .3)
ax.scatter(iter_x_nr,iter_y_nr,color = 'darkblue', marker = 'o',  label = 'Newton')
ax.quiver(iter_x_nr[:-1], iter_y_nr[:-1], anglesx_nr, anglesy_nr, scale_units = 'xy', angles = 'xy', scale = 1, color = 'darkblue', alpha = .3)
ax.legend()

ax.set_title('Comparing Newton and Gradient descent')

plt.show()