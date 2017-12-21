#!/usr/bin/env bash

wget -r -np -nH --user gec --password gecilluminadata --cut-dirs=100 -A .fastq.gz $1