import numpy as np
import torch
import gpytorch

def square(x):
    return x[0]**2 + x[1]

params = torch.tensor([1.,1.])
params.requires_grad_()

print(square(params))
# Do gradient descent
n_optim_steps = int(1e2)
optimizer = torch.optim.Adam([params], 1e-1)

for ii in range(n_optim_steps):
    optimizer.zero_grad()
    loss = square(params)

    print('Step # {}, loss: {}, parameter1: {}, parameter2: {}'.format(ii, loss.item(), params[0].item(), params[1].item()))
    loss.backward()
    # Access gradient if necessary
    grad = params.grad.data
    optimizer.step()
