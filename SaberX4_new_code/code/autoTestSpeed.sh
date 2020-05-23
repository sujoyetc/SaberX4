#!/bin/sh
make clean
make
cd saberX4Test
./test_kex
cd ..
