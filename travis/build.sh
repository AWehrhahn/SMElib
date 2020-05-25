#!/bin/bash

# Install dependencies
yum install -y gcc-gfortran cmake autoconf automake libtool

./bootstrap
./configure --prefix=$PWD

make install

# Only test with Python3
for PYBIN in /opt/python/cp3*/bin; do
    echo "${PYBIN}"
    "${PYBIN}/pip" install -r ./test/requirements.txt
    "${PYBIN}/pytest"
done

ls /io/lib
ls /io/share/smelib