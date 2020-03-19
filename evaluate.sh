#!/bin/bash

g++ main.cpp -std=c++17 -lemon -O3 -fopenmp

g_phi="0.05"
h_phi="0.4"
h_rat="0.02"
files="$@"

mkdir -p result

for f in $files; do
  time timeout 30m ./a.out --G_phi="$g_phi" --H_phi="$h_phi" --vol=1 --h_ratio="$h_rat" -f "$f" 2>&1| tee -a result/$(basename $f).out ; sleep 1s;
done
