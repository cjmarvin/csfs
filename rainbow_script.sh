#!/usr/bin/env sh
#
# This is a self-running script to produce a colorized echellogram.
#
# Two shell lines must be given:
# 1) Run CSFS to produce a spectrum to be colorized
# 2) Run CSFS to produce a wavemap

python csfs.py V 2 -mc -a P -b X --out --path --tell --noise --blaze -nr 1e10
python csfs.py V 2 -mc -a W -b W --out viswavemap -ns 1000