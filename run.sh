#!/usr/bin/sh

echo "Steps, x0, x1, x2, x3, x4, x5" > results.csv
for i in {19..19}
do
    ./sampling.exe $i | tee -a results.csv
done