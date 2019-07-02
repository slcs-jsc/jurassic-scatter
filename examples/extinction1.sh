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
$src/climatology extinction1.ctl - atm.tab || exit

info "Create geometry..."
$src/limb extinction1.ctl 800 5 15 1 obs.tab || exit

info "Call forward model..."
$src/formod extinction1.ctl obs.tab atm.tab rad_ext1.tab AEROFILE aero0.tab|| exit

# info "Compare files..."
compare rad_ext1.tab
# compare atm.tab
# compare obs.tab
# compare rad.tab.CO2
# compare rad.tab.EXTINCT
# compare rad.tab.H2O
# compare rad.tab.O3
