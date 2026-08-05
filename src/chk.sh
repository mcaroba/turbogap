#!/usr/bin/bash

file=$1

grep malloc $file | sed 's/!.*//g' |  cut -d '(' -f 2  | cut -d ',' -f 1 | sed '/^[[:space:]]*$/d' | sort > t1
grep gpu_free $file | sed 's/!.*//g' |  cut -d '(' -f 2  | cut -d ',' -f 1 | sed '/^[[:space:]]*$/d' |sort > t2

# grep malloc $file  | sort > t1
# grep gpu_free $file  |sort > t2

paste t1 t2 > tmp
rm t1 t2

