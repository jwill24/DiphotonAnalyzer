#!/bin/sh
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input files> <output file> [extra flags]"
    exit
fi
hadd $3 $2 `xrdfs root://eoscms.cern.ch ls -u $1 | grep '\.root'`
