#!/bin/bash
source ~/.bashrc
for fil in `ls *.f *.dek` Makefile; do echo $fil; mydiff $fil $1/$fil; done
