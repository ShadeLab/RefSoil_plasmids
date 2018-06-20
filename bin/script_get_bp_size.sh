#!/bin/bash -login

DATA=$1

head -n1 -q *.gbk | awk -v OFS="\t" '$1=$1' > ${DATA}_size_refsoil_tab.txt
