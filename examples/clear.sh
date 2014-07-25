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
$src/climatology clear-air.ctl - atm.tab CLIMZONE pwin|| exit

info "Create observation geometry..."
$src/limb clear-air.ctl 800 5 15 1 obs.tab || exit

info "Call forward model..."
$src/formod clear-air.ctl obs.tab atm.tab rad_clear.tab || exit

info "Compare files..."
compare atm.tab
compare obs.tab
compare rad_clear.tab
# compare rad.tab.CO2
# compare rad.tab.EXTINCT
# compare rad.tab.H2O
# compare rad.tab.O3
