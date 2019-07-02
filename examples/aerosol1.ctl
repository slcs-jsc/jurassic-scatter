# Emitters...
NG = 3
EMITTER[0] = CO2
EMITTER[1] = H2O
EMITTER[2] = O3

# Channels...
ND = 1
NU[0] = 792

# Tables...
TBLBASE = tab/boxcar

# Aerosol
SCA_N = 6
SCA_MULT = 1

MAX_QUEUE = 1e5

# Atmosphere/Climatology
DZ = 1
ZMAX = 80
ZMIN = 0 
CLIMPATH = ../clim/clim_remedios.tab

#diese Eingabe ist auch m√?glich, aber nicht wirklich n√?tig
#DIRLIST = dirlist-aero.asc
#AEROFILE = aero.tab
