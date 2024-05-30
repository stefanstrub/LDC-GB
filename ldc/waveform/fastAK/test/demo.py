import numpy as np
import matplotlib.pyplot as plt
import ldc.waveform.fastAK as fastAK
from ldc.common.tools import compute_tdi_snr
from ldc.common.series import TDI
from ldc.lisa.noise import get_noise_model

maxDur = 62914560.
timestep = 15.0
dtPhase = 2048.
df = 1.0/maxDur
Extended = 0

ak_fd = fastAK.pyFreqAK_RA(timestep, maxDur, dtPhase, 0)

params = {'MBHMass': 1338797.0898213838,
          'mu': 147.67680599627957, 'spin': 0.508191849430434,
          'lam': 0.5159259750304522, 'thS': 1.9963332443330797,
          'phS': 2.3328932152936743, 'thK': 1.3576019682706004,
          'phK': 4.464051088924649, 'DL': 26079065546.869873,
          'e0': 0.16954807413751763, 'nu0': 0.0014474837386423904,
          'phi0': 4.065749328184277, 
          'alph0': 5.26964680932448, 'gam0': 2.010025317378693,
          't0': 41097264.4945925}

Xf, Yf, Zf = ak_fd.get_fd_tdixyz(params)

plt.figure
plt.loglog(Xf.f, np.abs(Xf))
#plt.show()

noise = get_noise_model("SciRDv1", Xf.f)
noise.set_wdconfusion(2)

tdi = TDI([Xf, Yf, Zf], names=["X", "Y", "Z"])
tdi.XYZ2AET()
snr = compute_tdi_snr(tdi.ds, noise, AET=True)

print(snr)
print(np.sqrt(snr["A2"]))
#SNRs: 12.526737668034334 14.653158441218125 0.02961360655254664
