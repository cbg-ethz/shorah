#!/bin/sh -l

autoreconf -vif -I m4
./configure 
make -j1 distcheck