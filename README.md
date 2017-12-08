# simuHybrid
Scripts (.py scripts and wrappers) for simulating hybrid swarms using the Python simuPOP environment. The goal of these scripts is to simulate hybridization and recombination within populations to track the fate of alleles at two types of loci: pairs of loci that interact through epistasis and three loci taht additively affect fitness in the environment. 

As indicated in the wrapper scripts, parameters that can be varied include 1) s_dmis: the strength of selection acting on the epistatic loci, 2) s_adds: the strength of selection acting on 'adaptive' loci (it's a bad term for those loci, but it is what it is), 3) rs: recombination rates between adjacent loci (actual rate is 2x the number specified) 4) ns: population sizes within the hybrid zone and the two parental demes 5) ms: probabilities of migrating into the hybrid zone (mean number of immigrants from each parental population is N * m).

Two aspects of the model are easily flexible, but currently hard-coded into the scripts. First, selection acting on individuals as a function of their genotype at epistatic and adaptive loci can take one of three forms: (1) "sel_1" - selection against mixed ancestry within epistatic pairs and selection against '0' alleles at adaptive loci; ie. 'directional selection' model. (2) "sel_2" - selection against mixed ancestry within epistatic pairs and selection favoring mixed ancestry at adaptive loci; ie. 'selection-for-admixture' model. (3) "sel_3" - selection against mixed ancestry within epistatic pairs and selection against admixture at adaptive loci; ie. 'disruptive selection' model. To change between models of selection, one must manually specify the desired "sel" function on line 234 of the script (where the PySelector is called). Second, the starting genetic structure of the hybrid swarm can be either all F1 hybrids or a mixture of randomly mating parental genotypes. To switch between initial conditions one has to manually comment in or out the appropriate lines between lines 201-206.

Each "simuHybrid_hybrids_w_geneflow*" script then runs simulations of hybrid swarms, records allele frequencies at loci that affect the fitness of an individual and are arranged in one of three different genetic architectures: "dispersed", "interspersed", or "modular".

Required Python libraries: simuPOP, numpy, scipy, pprint, subprocess (call), sys, time, glob, fileinput, re, math, os

<b>To Run:</b>
<p>1) download master branch and place all files into a single directory.

2) open a terminal window and cd into the directory containing the scripts.

3) make sure .py scripts are executable (e.g. chmod u+x <file_name.py>).

4) run wrapper scripts (e.g. ./submit_hybrids_w_geneflow_disperseArchi.py)</p>

  - there are three wrappers, one for each of three types of genetic architecture: "dispersed", "interspersed", and "modular".
  - running a given wrapper will call the corresponding call_*.sh script and submit jobs on a cluster running slurm. You will of course want to tweak these scripts to work in the specific environment you're operating in.
 
<b>Output:</b>
<p>.txt files that contain allele frequencies at simulated loci, recorded every 10 generations for a total of 1000 generations. The ordering of allele frequencies are: rows are allele frequencies of parent deme 1, hybrid deme, and parent deme 2. Columns are first epistatic pair, second epistatic pair, followed by the three adaptive loci.</p>

  - The header of the output files specifies the version of simupop ran, and the parameters used are given on the line starting with "argument successfully loaded:". The following blocks are allele frequencies for each deme (rows) and locus (columns) every 10 generations.
  - Each replicate simulation is given in sequential order. So after allele frequencies are reported for the 1000th generation of given iteration, the following lines will be allele frequencies at the end of generation 0 for the next iteration under the given set of parameters.
  - A seperate .txt file in generated for every different combination of parameter values hybrid swarms are simulated under.
  - The parameter values are specified in the file name. One should check that the parameter values reported in the file names match those specified in the wrappers.
