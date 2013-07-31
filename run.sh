#!/usr/bin/sh

echo "Steps, x0, x1, x2, x3, x4, x5" > results.csv
for i in {0..14}
do
    ./sampling.exe $i | tee -a results.csv
done