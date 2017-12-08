#!/usr/bin/python/2.7.12

'''

Submits jobs through 'call_simuHyb_w_geneflow_disperseArchi.sh' for all parameter combinations.

usage:
python submit_hybrids_w_geneflow_disperseArchi.py

'''


import sys, os

s_dmis = [0.000, 0.001, 0.01, 0.05, 0.1]
s_adds = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
h = .5
rs = [0.5, 0.1, 0.05] 
ns = [50, 100, 1000]
ms = [0.0001, 0.001, 0.01]


for n in ns:
    for r in rs:
        for m in ms:
            for s_dmi in s_dmis:
                for s_add in s_adds:
                    
                    os.system( 'sbatch call_hybrids_w_geneflow_disperseArchi.sh {0} {1} {2} {3} {4} {5}'.format(  s_dmi, s_add, h, r, n, m ) )
    
