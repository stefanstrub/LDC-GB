import numpy as np

opticalnoise = 5.00e-12 #
readoutnoise = 6.35e-12 #
accnoise =  2.4e-15 # m/s^2/sqrt(Hz)

#seed = ...
#np.random.seed(seed)

#(1+a)^2 Shf = (1+b)^2 Sro + (1+b)^2 Sop
#(1+a)^2 Shf = (1+b)^2 ( Sro + Sop )
# 20% pertubation of Shf (a=0.2) ~ 20% pertubation of Sro and Sop (b=0.2)

for v,p in zip([opticalnoise, readoutnoise, accnoise],
               [0.2, 0.2, 0.2]):

    v += v*np.random.randn(1)*p

    print(v)

    
