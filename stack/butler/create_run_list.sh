#!/bin/bash

#ls -1c /sps/lsst/groups/FocalPlane/SLAC/run5/

input="run_list.txt"
while read -r line
do
  echo -n "'$line' "
done < "$input"
