#!/bin/bash

cat *.cpp > junk
cat *.c >> junk
cat junk | listErrors $1
rm junk

# $1 = common or 2D or 3D

if [[ $1 = "common" ]]; then
  echo "numbering errors:"
  grep -r "ERROR:"
  grep -r "ERROR0:"
  grep -r "ERROR2:"
  grep -r "ERROR3:"
  grep -r "ERROR4:"
  grep -r "ERROR5:"
  grep -r "ERROR6:"
  grep -r "ERROR7:"
  grep -r "ERROR8:"
  grep -r "ERROR9:"
fi

if [[ $1 = "2D" ]]; then
  echo "numbering errors:"
  grep -r "ERROR:"
  grep -r "ERROR0:"
  grep -r "ERROR1:"
  grep -r "ERROR3:"
  grep -r "ERROR4:"
  grep -r "ERROR5:"
  grep -r "ERROR6:"
  grep -r "ERROR7:"
  grep -r "ERROR8:"
  grep -r "ERROR9:"
fi 

if [[ $1 = "3D" ]]; then
  echo "numbering errors:"
  grep -r "ERROR:"
  grep -r "ERROR0:"
  grep -r "ERROR1:"
  grep -r "ERROR2:"
  grep -r "ERROR4:"
  grep -r "ERROR5:"
  grep -r "ERROR6:"
  grep -r "ERROR7:"
  grep -r "ERROR8:"
  grep -r "ERROR9:"
fi

