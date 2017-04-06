#!/usr/bin/env bash

make all

CVOUTDIR="./arp100120"

rm -rf ${CVOUTDIR}

mkdir ${CVOUTDIR}

cd ${CVOUTDIR}

../p1 500 500

../p2 0 0 1 1 ../test-images/test0.ppm 20.png
../p2 0.5 0 1 0.5 ../test-images/test1.ppm 21.png
../p2 0 0 1 1 ../test-images/test2.ppm 22.png
../p2 0 0 0.5 0.5 ../test-images/test3.ppm 23.png
../p2 0 0 1 1 ../test-images/test4.ppm 24.png

../p3 0 0 1 1 ../test-images/test0.ppm 30.png
../p3 0.5 0 1 0.5 ../test-images/test1.ppm 31.png
../p3 0 0 1 1 ../test-images/test2.ppm 32.png
../p3 0 0 0.5 0.5 ../test-images/test3.ppm 33.png
../p3 0 0 1 1 ../test-images/test4.ppm 34.png

../p4 0 0 1 1 ../test-images/test0.ppm 40.png
../p4 0.5 0 1 0.5 ../test-images/test1.ppm 41.png
../p4 0 0 1 1 ../test-images/test2.ppm 42.png
../p4 0 0 0.5 0.5 ../test-images/test3.ppm 43.png
../p4 0 0 1 1 ../test-images/test4.ppm 44.png

echo "These images were generated with the following git sha1:" > sha1.txt
git rev-parse HEAD >> sha1.txt
