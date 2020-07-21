# -*- coding: utf-8 -*-
# pylint: disable=C0302

"""
Define the configuration for the simulation.

The parameters defined here are not paremeters on nodes or graphs, because
they change the fundamental structure of the simulation graph. The graph
must therefore be compiled by LISANode again if they are changed.

Authors:
    Jean-Baptiste Bayle <bayle@apc.in2p3.fr>

"""

## LISA Simulation Options

# Measurement sampling frequency, in Hz
LISA_MEASUREMENT_FS = 3.0

# Measurement downsampling factor
# This is the ratio LISA_PHYSICS_FS / LISA_MEASUREMENT_FS)
LISA_MEASUREMENT_DOWNSAMPLING = 1

# Sampling frequency at which physics is simulated, in Hz
# Deduced from LISA_MEASUREMENT_FS and LISA_MEASUREMENT_DOWNSAMPLING
LISA_PHYSICS_FS = LISA_MEASUREMENT_FS * LISA_MEASUREMENT_DOWNSAMPLING

# Upsampling factor used for orbit simulation
# This is the ratio LISA_PHYSICS_FS / LISA_ORBIT_FS)
LISA_ORBIT_UPSAMPLING = 1

# TDI sampling frequency, in s
LISA_TDI_FS = 5.0

# Sampling frequency used to compute the orbits, in Hz
# Deduced from LISA_ORBIT_UPSAMPLING and LISA_PHYSICS_FS
LISA_ORBIT_FS = LISA_PHYSICS_FS / LISA_ORBIT_UPSAMPLING

# Test mass acceleration noise shape
LISA_ACC_NOISE_A_LEVEL = 3E-15 #2.4E-15 #m/s^2/sqrt(Hz) 
LISA_ACC_NOISE_F_KNEE = 0.4E-3 #Hz

# Type of orbit used to compute inter-spacecraft travel times
# For Keplerian or from-file orbits, we compute travel times using relativistic
# corrections, which we bypass for polynomial travel times
# Must be one of the following:
#   * 'keplerian'         for analytic Keplerian orbits (see node `KeplerianOrbits`),
#   * 'from_file'         to read spacecraft positions and velocities from an external file
#   * 'polynomial_tt'     to bypass orbit and use a polynomial function of given order (0, 1 or 2,
#                         can be set with parameter 'polynomial_tt_order') as the travel times;
#                         polynomial coeffficients can be set using parameters Li, dLi, ddLi
#   * 'ldc'               to use LDC orbits and travel-time computation - LISANode is required to
#                         be run inside the LDC pipeline to access necessary shared dependencies
#
# Please note that this does not set the TDI generation
LISA_ORBIT_TYPE = 'keplerian'

# Path to GW strain file in HDF5 format, or None
# Projected strains must be given at receiver time, in a single dataset with one column per link
LISA_GW_INPUT = None

# Central frequency of the lasers, in MHz
# All laser offsets are relative offsets to this value
LISA_LASER_CENTRAL_FREQUENCY = 2.816E8

# Whether to apply doppler shifts to small frequency fluctuations (laser, sideband) and large
# frequency offsets (laser_frequency, sideband_frequency) used for noise coupling in LISA
LISA_SIMULATE_DOPPLER = False

# Whether frequency fluctuations of interferometric measurements should be published as outputs
LISA_PUBLISH_MEASUREMENTS = True

# Whether large frequency offsets of interferometric measurements should be published as outputs
LISA_PUBLISH_BEATNOTE_FREQUENCIES = True


## TDI Offline Processing Options

# Sampling frequency of measurements used as inputs, in Hz
# This value must match LISA_MEASUREMENT_FS
TDI_MEASUREMENT_FS = 3.0

# Input file format for TDIFromFile, must be 'text' or 'hdf5'
TDI_INPUT_FORMAT = 'hdf5'

# Generation of TDI algorithm, must be 1.5 or 2
# This only determines which TDI combination is computed
TDI_GENERATION = 1.5

# Whether estimated Dopplers computed in TDI should be published as outputs
TDI_PUBLISH_DOPPLER_ESTIMATES = True

# Whether intermediary variables should be published as outputs
TDI_PUBLISH_INTERVAR = True

# Whether Michelson X, Y and Z variables should be published as outputs
TDI_PUBLISH_MICHELSON = True

# Whether quasi-orthogonal A, E and T variables should be published as outputs
TDI_PUBLISH_ORTHVAR = True

# Whether Sagnac α, β, and ɣ variables should be publihsed as outputs
TDI_PUBLISH_SAGNAC = False

# Whether clock correction should be carried out on Michelson variables
TDI_CLOCKCOR_MICHELSON = True

# Whether clock correction should be carried out on quasi-orthogonal variables
TDI_CLOCKCOR_ORTHVAR = True

# Whether clock correction should be carried out on Sagnac variables
TDI_CLOCKCOR_SAGNAC = False

# Whether to use exact clock-noise corrections [1]
# [1] Clock-jitter reduction in LISA time-delay interferometry combinations, arXiv:2005.02430
TDI_USE_UPDATED_CLOCKCOR = True

# Whether to compensate for Doppler shifts
TDI_USE_DOPPLER_CORRECTION = False
