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
$src/limb clear-air.ctl 800 8 8 1 obs_raytrace.tab || exit

info "Call raytrace module..."
$src/raytrace clear-air.ctl obs_raytrace.tab atm.tab || exit

info "Compare files..."
compare atm.tab
compare obs_raytrace.tab
compare los.0
