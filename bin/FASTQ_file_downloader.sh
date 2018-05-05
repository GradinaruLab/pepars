#!/usr/bin/env bash

wget -r -np -nH --user gec --password "$2" --cut-dirs=100 -A .fastq.gz $1
