#! /bin/bash

# ----------------------------------------------------------------------

function compare {
    
    # Write info...
    echo -n "Compare $1: "
    
    # Compare...
    diff -q $1 org/$1 && echo "OK" || exit
    
    # Remove...
    rm $1
}

# ----------------------------------------------------------------------

function info {
    
    # Write info...
    echo
    echo "=================================================="
    echo $1
    echo "=================================================="
    echo
}

# ----------------------------------------------------------------------

# Setup...
src=../src

info "Create atmosphere..."
$src/climatology aerosol1.ctl - atm.tab || exit

info "Create geometry..."
$src/limb aerosol1.ctl 800 5 15 1 obs.tab || exit

info "Call forward model..."
$src/formod aerosol1.ctl obs.tab atm.tab rad_aero1.tab AEROFILE aero1.tab|| exit

# info "Compare files..."
compare rad_aero1.tab

