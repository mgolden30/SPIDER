

import numpy as np
from scipy.io import savemat

density = np.load('coarse_grained_density.npy')
density_r = np.load('rg_coarse_grained_density.npy')


print(density_r)

savemat("data.mat", {"rho": density, "rho_rg": density_r})