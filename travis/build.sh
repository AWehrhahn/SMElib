#!/bin/bash

# Install dependencies
yum install -y gcc f2c cmake autoconf automake libtool

f2c -w -a -C++ -Nn1604 -Nq1200 -dsrc/eos/ src/eos/*.f
f2c -w -a -c++ -dsrc/sme/ src/sme/*.f
rm makefile.am
mv makefile_f2c.am makefile.am

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
ls /io/share/libsme