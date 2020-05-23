#!/bin/sh
make clean
make
cd saberX4Test
rm *.txt
./PQCgenKAT_kem

