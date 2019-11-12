#! /bin/sh
# script to run code with some 
# small tsp library instances (<= 100)

rm -f results.csv
for inst in *.dist;
do
    ./tsp-subtour $inst
done
