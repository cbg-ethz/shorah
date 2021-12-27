#!/bin/bash

g++ -Wall -Wextra b2w.cpp -o out -lhts
g++ -g -Wall -Wextra -Iinclude/ src/dpm_sampler.cpp -lhts -DNDEBUG