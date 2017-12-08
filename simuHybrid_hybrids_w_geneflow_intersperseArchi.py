#!/usr/bin/python
#./working_test.py 0.1 0.1 0.5 0.1 100 0.01
'''

Forward time simulations of hybrid populations with migration with parental populations with simuPOP

#
#
# four loci as part of two pairs of epistatic "DMI" loci.
# three loci affecting a trait that additively affects fitness. 
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
import random

if len(sys.argv[1:]) == 6:
    vals = sys.argv[1:]
    print "\narguments successfully loaded: s.dmi =", vals[0], "s.add =", vals[1], "h =", vals[2], "r =", vals[3], "N=", vals[4], "m=", vals[5]
else:
    print "\nUsage: <simuHybrid_schumer_model_modularArchi.py> s_dmi(float) s_add(float) h(float) r(float) n(int)" 


# ---- user defined command line arguments: selection, dominance, recombination rate, and popn size  ----

s_dmi    = float(vals[0]) # selection on epistatic loci
s_add    = float(vals[1]) # selection on loci with additive effects on fitness
h        = float(vals[2]) # dominance coefficient at epistatic loci
r        = float(vals[3]) # recombination rate
popsize  = int(vals[4])   # population size
m        = float(vals[5]) # migration rate


# ---- define the position of loci involved in a given interaction and subject to a given type of selection ---- #
# location of epistatic loci:
dmi1 = 0  
dmi2 = 6
dmi3 = 12
dmi4 = 18

# location of loci with additive fitness effects:
ad1 = 3
ad2 = 9
ad3 = 15




##########################################################
## set up functions to define fitness based on genotype ##
##########################################################
def sel_1(geno):
	# selection is against mixed ancestry at 'epistatic' loci and against '0' alleles at 'adaptive' loci #
	
	# 1) 'good' dmi genotypes:
    if sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2, dmi3, dmi4]] ) == 8 or \
    	sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2, dmi3, dmi4]] ) == 0 or \
       	sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 4 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 0 or \
       	sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 0 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 4:
       	return 1 - ((6 - sum( [sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]] )) * s_add) 
    
    # 2) genotypes w/ 1 mismatched allele:
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 4 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2 or \
    	sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 0 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2 or \
    	sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 4 or \
    	sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 0:
    	return 1- (s_dmi) - ((6 - sum( [sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]] )) * s_add) 
    
    # 3) both pairs carrying 2 'bad' alleles of 4 possible alleles:
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 2:
    	return 1 - (2*s_dmi) - (2*s_dmi) - ((6 - sum( [sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]] )) * s_add)
    
    # 4) the scenario of genotypes where both dmi pairs have 1 missmatched allele.
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2: 
        return 1 - (s_dmi) - (s_dmi) - ((6 - sum( [sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]] )) * s_add)
    
    # 5) the scenario w/ 1 pair carrying 2 missmatched alleles and the other pair carrying one missmatched allele:
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2 or \
    	sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 2:
    	return 1 - (2*s_dmi) - (s_dmi) - ((6 - sum( [sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]] )) * s_add)
    
    # 6) one pair w/ 'good' interactions, the other carrying two 'bad' allele:
    else:
    	return 1 - (2*s_dmi) - ((6 - sum( [sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]] )) * s_add)




# 2) selection for admixture - e.g. transgressive segregation
def transgressive_sel(geno, single_geno):
        if sum([sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]])%6 == 0:
                return 2
        if sum([sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]]) == 1 or sum([sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]]) == 5:
                return 2 - sum(single_geno)%2
        else:
                return sum(single_geno)%2

def sel_2(geno):
        # selection is against mixed ancestry at 'epistatic' loci and admixed haplotypes are favored at 'adaptive' loci #
    
        # 1) 'good' dmi genotypes:
    if sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2, dmi3, dmi4]] ) == 8 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2, dmi3, dmi4]] ) == 0 or \
        sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 4 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 0 or \
        sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 0 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 4:
        return 1 - (sum([transgressive_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 2) genotypes w/ 1 mismatched allele:
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 4 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 0 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 4 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 0:
        return 1- (s_dmi) - (sum([transgressive_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 3) both pairs carrying 2 'bad' alleles of 4 possible alleles:
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 2:
        return 1 - (2*s_dmi) - (2*s_dmi) - (sum([transgressive_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 4) the scenario of genotypes where both dmi pairs have 1 missmatched allele.
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2:
        return 1 - (s_dmi) - (s_dmi) - (sum([transgressive_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 5) the scenario w/ 1 pair carrying 2 missmatched alleles and the other pair carrying one missmatched allele:
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 2:
        return 1 - (2*s_dmi) - (s_dmi) - (sum([transgressive_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 6) one pair w/ 'good' interactions, the other carrying two 'bad' allele:
    else:
        return 1 - (2*s_dmi) - (sum([transgressive_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)




# 3) disruptive selection favoring parental genotypes:
def disrupt_sel(geno, single_geno):
        if sum([sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]])%6 == 0:
                return 0
        elif sum([sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]]) == 3 and sum(single_geno)%2 == 0:
                return 2
        elif sum([sum(geno[(ad*2):(ad*2)+2]) for ad in [ad1, ad2, ad3]]) < 3:
                return sum(single_geno)
        else:
                return 2-sum(single_geno)

def sel_3(geno):
        # selection is against mixed ancestry at 'epistatic' loci and parental haplotypes are favored at 'adaptive' loci #

        # 1) 'good' dmi genotypes:
    if sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2, dmi3, dmi4]] ) == 8 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2, dmi3, dmi4]] ) == 0 or \
        sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 4 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 0 or \
        sum(geno[dmi1*2:(dmi1*2)+2]) + sum(geno[dmi2*2:(dmi2*2)+2]) == 0 and sum(geno[dmi3*2:(dmi3*2)+2]) + sum(geno[dmi4*2:(dmi4*2)+2]) == 4:
        return 1 - (sum([disrupt_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 2) genotypes w/ 1 mismatched allele:
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 4 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 0 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 4 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 0:
        return 1- (s_dmi) - (sum([disrupt_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 3) both pairs carrying 2 'bad' alleles of 4 possible alleles:
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 2:
        return 1 - (2*s_dmi) - (2*s_dmi) - (sum([disrupt_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 4) the scenario of genotypes where both dmi pairs have 1 missmatched allele.
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2:
        return 1 - (s_dmi) - (s_dmi) - (sum([disrupt_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 5) the scenario w/ 1 pair carrying 2 missmatched alleles and the other pair carrying one missmatched allele:
    elif sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) == 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) != 2 or \
        sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi1, dmi2]] ) != 2 and sum( [sum(geno[(dmi*2):(dmi*2)+2]) for dmi in [dmi3, dmi4]] ) == 2:
        return 1 - (2*s_dmi) - (s_dmi) - (sum([disrupt_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)

    # 6) one pair w/ 'good' interactions, the other carrying two 'bad' allele:
    else:
        return 1 - (2*s_dmi) - (sum([disrupt_sel(geno, geno[pos*2:(pos*2)+2]) for pos in [ad1,ad2,ad3]]) * s_add)




###################
# run simulations #
###################

for rep in range(1,501):
    pop = sim.Population([popsize]*3, ploidy = 2, loci=[20], lociPos=[1,1.1,2,3,3.1,4,5,5.1,6,7,7.1,8,9,9.1,10,11,11.1,12,13,13.1], infoFields=['fitness','migrate_to'])
    
    # initiate genotypes in subpopulation 0 as fixed for '0' ancestry:
    sim.initGenotype(pop, genotype=[0]*20, subPops=[0])
    
    # genotypes at subpopulation 1:
    #sim.initGenotype(pop, genotype=[0]*20+[1]*20, subPops=[1])    # all F1 hybrids
    
    sim.initGenotype(pop, genotype=[0]*20, subPops=[1])            # or a mixture of 
    for ind in random.sample(range(popsize,popsize*2), popsize/2): # both parental species.
        pop.individual(ind).setGenotype([1]*20)
    
    # initiate genotypes in subpopulation 2 as fixed for '1' ancestry
    sim.initGenotype(pop, genotype=[1]*20, subPops=[2])
    
    def printAlleleFreq(pop):
        'Print allele frequencies of all loci and populations'
        
        sim.stat(pop, alleleFreq=[dmi1, dmi2, dmi3, dmi4, ad1, ad2, ad3], vars=['alleleFreq_sp'])
        
        print 'Allele frequencies at generation', pop.dvars().gen
        
        for p in range(3):
            for l in [dmi1, dmi2, dmi3, dmi4, ad1, ad2, ad3]:
                if l == ad3:
                    print '%.2f\n' % pop.dvars(p).alleleFreq[l][1],
                else:
                    print '%.2f' % pop.dvars(p).alleleFreq[l][1],
        
        return True
    
    pop.evolve(
        initOps = [
            sim.InitSex(),
            #sim.Stat(popSize=True),
            #sim.PyEval(r'"%d %s " % (gen, subPopSize)'),
        ],
        
        preOps=sim.PySelector(loci=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19], func=sel_1),  # !! change for different types of selection !! #
        
        matingScheme = sim.RandomMating(ops=sim.Recombinator(intensity=r), subPopSize=[popsize,popsize,popsize]),
        
        postOps = [
            sim.Migrator(rate=
                [
                [0.0, m, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, m, 0.0]
                ]
            ),
            #sim.Stat(popSize=True),
            #sim.PyEval(r'"%d %s " % (gen, subPopSize)'),
            sim.PyOperator(printAlleleFreq, step=10),
        ],
        gen = 1001
    )
