from few.utils.utility import check_for_file_download
import os
import few

fp = "AmplitudeVectorNorm.dat"
fp = "FluxNewMinusPNScaled_fixed_y_order.dat"
fp = "SchwarzschildEccentricInput.hdf5"
fp = "Teuk_amps_a0.0_lmax_10_nmax_30_new.h5"

dir_path = os.path.dirname(os.path.realpath(few.__file__))
few_dir = dir_path + '/../'

for fp in ["AmplitudeVectorNorm.dat", "FluxNewMinusPNScaled_fixed_y_order.dat",
           "SchwarzschildEccentricInput.hdf5",#"Teuk_amps_a0.0_lmax_10_nmax_30_new.h5"
           ]:
    check_for_file_download(fp, few_dir)
