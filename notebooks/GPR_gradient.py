import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF,ConstantKernel
from scipy.special import erf

#Gives a fitted Gaussian Process object that can then be used for predictions.
#The Input is of the Form x.shape = (n), y.shape = (n,t) where both x and y
#are np.ndarrays.
#The normalisation has to be set to False for now since it didn't work with
#my current version of sklearn. Could be added in customary by normalizing the
#input data and denormalizing the output directly.
#The Kernel types (not their parameters though) have to stay this way since the derivates
#and antiderivates are computed for this setup. Should no constant kernel be 
#desired its parameters can be set to constant_value = 1.0 and 
#constant_value_bounds = 'fixed'.
#All other values, as n_restarts, the RBF kernel and Constant kernel parameters
#have to be selected according to the input data.

class GPR:
    def __init__(self,x,y):
        normalize = False #hardcoded, don't change.
        n_restarts = 2

        # k1 = ConstantKernel(constant_value=1.0,constant_value_bounds=(1e-5,1e5))
        k2 = RBF(length_scale=0.1,length_scale_bounds=(1e-5,1e5))

        self.gpr = GaussianProcessRegressor(k2,
                                           n_restarts_optimizer=n_restarts,
                                           normalize_y=normalize).fit(x.reshape(-1,1),y)

    def predict(self,x,k=0):
        #x of shape (m)
        
        #returns the gp predictions where f is the true function and
        #df, ddf, If, IIf are its first and second derivate respectively antiderivates
        #the outputs are the predictions f_p,df_p,ddf_p,If_p,IIf_p where
        #f(x) = f_p(x), df(x) = df_p(x), ddf(x) = ddf_p(x), If(x) = If_p(x) + C1, 
        #IIf(x) = IIf_p(x) + C1*x + C2 with some constants C1,C2
        #set k = 0 for the normal prediction, K = 1,2 for the first or second derivates
        #and k = -1,-2 for the first or second antiderivates
    
        x = x.reshape(-1,1)
    
        X = x - self.gpr.X_train_.reshape(1,-1)
        # c = self.gpr.kernel_.k1.constant_value
        l = self.gpr.kernel_.length_scale
        A = self.gpr.alpha_

        f = np.exp(-(X)**2 / (2*l**2))
        df = (f * (-X / l ** 2))

        if k == 0: 
            return  f @ A
        elif k == 1: 
            return  df @ A
        else:
            raise Exception('Unknown parameter k: {}'.format(k))

    def predict2(self, x):
        # gets 'l' used in denominator of expected value of gradient for RBF kernel 
        k2_l = self.gpr.kernel_.length_scale

        # not necessary to do predict, but now y_pred has correct shape
        y_pred, sigma = self.gpr.predict(np.asanyarray([x]).reshape(-1,1),  return_std=True)

        # allocate array to store gradient
        y_pred_grad = 0.0*y_pred

        # set of points where gradient is to be queried
        # x = np.atleast_2d(np.linspace(-5, 0.8, 1000)).T
        
        X = self.gpr.X_train_
        # loop over each point that a gradient is needed
        for key, x_star in enumerate(x):
            # eval_gradient can't be true when eval site doesn't match X
            # this gives standard RBF kernel evaluations
            k_val= self.gpr.kernel_(X, np.atleast_2d(x_star), eval_gradient=False).ravel()

            # x_i - x_star / l^2
            x_diff_over_l_sq = ((X-x_star)/np.power(k2_l,2)).ravel()

            # pair-wise multiply
            intermediate_result = np.multiply(k_val, x_diff_over_l_sq)

            # dot product intermediate_result with the alphas
            final_result = np.dot(intermediate_result, self.gpr.alpha_)

            # store gradient at this point
            y_pred_grad[key] = final_result
        return y_pred_grad

def function(x):
    return np.sin(x)

def dfunction(x):
    return np.cos(x)

x_data = np.random.rand(100)*np.pi
y_data = function(x_data)

train_test_ratio = 0.8
n_split = int(len(x_data)*train_test_ratio)
x_train = x_data[:n_split]
y_train = y_data[:n_split]
x_test = x_data[n_split:]
y_test = y_data[n_split:]

GPR1 = GPR(x_train, y_train)
y_predicted = GPR1.predict(x_test)
RMSE = np.sqrt(np.mean((y_predicted - y_test)**2))

yd_predicted = GPR1.predict(x_test, k=1)
yd_predicted2 = GPR1.predict2(x_test)
RMSE_grad = np.sqrt(np.mean((yd_predicted - dfunction(x_test))**2))
RMSE_grad2 = np.sqrt(np.mean((yd_predicted2 - dfunction(x_test))**2))

print('end')