#! /bin/bash

# ----------------------------------------------------------------------
function infomsg {
# ----------------------------------------------------------------------
    
    echo
    echo "============================================================"
    echo "Compile: $1"
    echo "============================================================"
    echo
}

# ----------------------------------------------------------------------
# Main...
# ----------------------------------------------------------------------

# Check arguments...
if [ $# -ne 1 ] ; then
    echo "usage: $0 <target-dir>"
    exit
fi

# Set number of cores...
THREADS=$(cat /proc/cpuinfo | grep processor | wc -l)

# Get absolute path...
target=$(mkdir -p $1 && cd $1 && pwd)

# Prepare directories...
mkdir -p $target/src $target/bin $target/lib $target/man/man1 \
    && cp *tar.bz2 $target/src \
    && cd $target/src \
    && for f in $(ls *tar.bz2) ; do tar xvjf $f ; done \
    || exit

# GSL...
dir=gsl-1.15
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make -j$THREADS && make check && make install && make clean \
    || exit

# netCDF...
dir=netcdf-4.1.3
infomsg $dir
cd $target/src/$dir \
    && ./configure --prefix=$target --disable-netcdf-4 --disable-dap --disable-fortran --disable-cxx \
    && make -j$THREADS && make check && make install && make clean \
    || exit
