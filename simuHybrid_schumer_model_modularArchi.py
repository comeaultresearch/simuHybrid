#!/usr/bin/python
#./simuHybrid_schumer_model_modularArchi.py 0.1 0.1 0.5 0.1 50
'''

Forward time simulations of hybrid swarms with simuPOP

#
#
# four loci as part of two pairs of epistatic "DMI" loci.
# three loci affecting a trait that additively affects fitness. 
# 
# 
#  
 
'''

# ---- overhead and usage information ----
import simuPOP as sim
import numpy as np
import scipy as sp
from pprint import pprint
from subprocess import call
import sys, time, glob, fileinput, re, math, os


if len(sys.argv[1:]) == 5:
    vals = sys.argv[1:]
    print "\narguments successfully loaded: s.dmi =", vals[0], "s.add =", vals[1], "h =", vals[2], "r =", vals[3], "N=", vals[4]
else:
    print "\nUsage: <simuHybrid_schumer_model_modularArchi.py> s_dmi(float) s_add(float) h(float) r(float) n(int)" 


# ---- user defined command line arguments: selection, dominance, recombination rate, and popn size  ----

s_dmi    = float(vals[0]) # selection on epistatic loci
s_add    = float(vals[1]) # selection on loci with additive effects on fitness
h        = float(vals[2]) # dominance coefficient at epistatic loci
r        = float(vals[3]) # recombination rate
popsize  = int(vals[4])        # population size


# ---- file handling / management  ----

# get the directory that the script has been run from:
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)

# specify path to output file:
my_dir_string = dname+'/simuHyb_modularArch_n'+str(popsize)+'_'+time.strftime("%d-%m-%Y")
call(["mkdir", my_dir_string])
out = open(my_dir_string+'/simuHyb_modularArch_'+time.strftime("%d-%m-%Y")+'_data_sDMI'+str(s_dmi)+'_sADD'+str(s_add)+'_r'+str(r)+'.txt', "w")

# write header line for modular genetic architecture:
out.write('gen\tdmi1\tdmi1_lnkd\tneu1\tdmi2\tdmi2_lnkd\tneu2\tdmi3\tdmi3_lnkd\tneu3\tdmi4\tdmi4_lnkd\tneu4\tadapt1\tadapt1_lnkd\tneu5\tadapt2\tadapt2_lnkd\tneu6\tadapt3\tadapt3_lnkd\n')
out.close()


## define the position of loci involved in a given interaction and subject to a given type of selection ##
# location of epistatic loci:
dmi1 = 0  
dmi2 = 3
dmi3 = 6
dmi4 = 9

# location of loci with additive fitness effects:
ad1 = 12
ad2 = 15
ad3 = 18


# ---- ---- 

# run simulations
for rep in range(1,6):
    pop = sim.Population(size=50, ploidy = 2, loci=[20], lociPos=[1,1.1,2,3,3.1,4,5,5.1,6,7,7.1,8,9,9.1,10,11,11.1,12,13,13.1], infoFields=['fitness'])
    sim.initGenotype(pop, genotype=[0]*20+[1]*20)  # all individuals are heterozygous across the genome
    
    #########################################################
    ## set up function to define fitness based on genotype ##
    #########################################################
    def sel(geno):

        # 'good' dmi genotypes:
        if sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 8 or \
           sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 0 or \
           sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 4 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 0 or \
           sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 0 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 4:
            return 1 - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add) 
        
        # second, the scenario of genotypes w/ 1 mismatched allele:
        elif sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 7 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 1:
            return 1- (h*s_dmi) - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add) 
       
        # third, the scenario of genotypes where both dmi pairs have 1 missmatched allele.
        elif sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 3 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 3 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 1 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 1: 
            return 1 - (2*h*s_dmi) - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add)
            
        # fourth, the scenario w/ 1 homozygous missmatched pair:
        elif sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 6 and sum(geno[dmi1*2:(dmi1*2)+2]) == 0 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 6 and sum(geno[dmi2*2:(dmi2*2)+2]) == 0 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 6 and sum(geno[dmi3*2:(dmi3*2)+2]) == 0 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 6 and sum(geno[dmi4*2:(dmi4*2)+2]) == 0 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 2 and sum(geno[dmi1*2:(dmi1*2)+2]) == 2 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 2 and sum(geno[dmi2*2:(dmi2*2)+2]) == 2 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 2 and sum(geno[dmi3*2:(dmi3*2)+2]) == 2 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 2 and sum(geno[dmi4*2:(dmi4*2)+2]) == 2:
            return 1-(2*s_dmi) - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add)
        
        # fifth, the scenario w/ 1 heterozygous pair:
        elif sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 6 and sum(geno[dmi1*2:(dmi1*2)+2]) == 1 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 6 and sum(geno[dmi3*2:(dmi3*2)+2]) == 1 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 2 and sum(geno[dmi1*2:(dmi1*2)+2]) == 1 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 2 and sum(geno[dmi3*2:(dmi3*2)+2]) == 1 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 4 and sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 3 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 4 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 3:
            return 1 - (2*h*s_dmi) - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add)

        # sixth, the scenario w/ 1 missmatched allele, but alternately fixed DMI pairs:
        elif sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 4 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 1 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 1 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 4 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 0 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 3 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 3 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 0:
            return 1 - (h*s_dmi) - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add)
        
        # seventh, the scenario with one pair carrying a single missmatched allele and the other homozygous missmatched e.g. 01/02 or 21/02
        elif sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 5 and sum(geno[dmi3*2:(dmi3*2)+2]) != 1 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 5 and sum(geno[dmi1*2:(dmi1*2)+2]) != 1 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 3 and sum(geno[dmi3*2:(dmi3*2)+2]) != 1 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 3 and sum(geno[dmi1*2:(dmi1*2)+2]) != 1:
            return 1 - (h*s_dmi) - (2*s_dmi) - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add)
  
        # eighth, the scenario w/ 3 heterozygous loci:
        elif sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 5 or \
             sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 3:       
           return 1 - (h*s_dmi) - (2*(h*s_dmi)) - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add)
   
        # ninth, heterozyous at all 4 dmi loci:
        elif sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) + sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 4 and sum(geno[dmi1*2:(dmi1*2)+2]) == 1 :
           return 1 - (4*(h*s_dmi)) - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add)
        
        # finally, homozyous missmatched at both DMI pairs:
        else:
           return 1 - (2*s_dmi) - (2*s_dmi) - ((6 - (sum(geno[ad1*2:(ad1*2)+2]) + sum(geno[ad2*2:(ad2*2)+2]) + sum(geno[ad3*2:(ad3*2)+2]))) * s_add)

    
    
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.PyOutput("starting 'pop.evolve'\n"),
        ],
        
        preOps=sim.PySelector(loci=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19], func=sel), 
        matingScheme=sim.RandomMating(ops=sim.Recombinator(intensity=r)),
        
        postOps=[
            sim.Stat(alleleFreq=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19], step=100),
            sim.PyEval(r"'%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' % "
                "(gen, alleleFreq[0][1], alleleFreq[1][1], alleleFreq[2][1], "
                "alleleFreq[3][1], alleleFreq[4][1], alleleFreq[5][1], "
                "alleleFreq[6][1], alleleFreq[7][1], alleleFreq[8][1], "
                "alleleFreq[9][1], alleleFreq[10][1], alleleFreq[11][1], "
                "alleleFreq[12][1], alleleFreq[13][1], alleleFreq[14][1], "
                "alleleFreq[15][1], alleleFreq[16][1], alleleFreq[17][1], "
                "alleleFreq[18][1], alleleFreq[19][1])", step=100,
                output = '>>>'+my_dir_string+'/simuHyb_modularArch_'+time.strftime("%d-%m-%Y")+'_data_sDMI'+str(s_dmi)+'_sADD'+str(s_add)+'_r'+str(r)+'.txt')
        ], 
        gen = 1001
    )
