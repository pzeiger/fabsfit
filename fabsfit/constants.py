import scipy.constants as cst


hbar = cst.physical_constants['reduced Planck constant times c in MeV fm'][0] * \
       1000 * 1e-5    # MeV -> keV and fm -> Å
me = cst.physical_constants['electron mass energy equivalent in MeV'][0] * \
       1000           # MeV -> keV

c = cst.physical_constants['speed of light in vacuum'][0]  # in m/s
