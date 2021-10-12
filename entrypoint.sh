#!/bin/sh -l

autoreconf -vif -I m4
./configure 
make install
make -j1 distcheck