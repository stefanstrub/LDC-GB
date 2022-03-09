#!/bin/sh
for i in `seq 0 13`
do
  echo "run batch: $i"
  python notebooks/LDC1-4_even3.py $i
done

