# simuHybrid
Scripts (.py scripts and .R wrappers) for simulating hybrid swarms using the Python simuPOP environment. Goal is to track allele frequency change within hybrid swarms at two type of loci: epistatic pairs and those additively affecting fitness in the environment. Hybrid swarms begin as a pool of F1s.

Each "simuHybrid_schumer_model_*" script runs simulations of hybrid swarms, tracking allele frequencies at loci that affect the fitness of an individual and are arranged in one of three different genetic architectures: "dispersed", "interspersed", or "modular".

Required Python libraries: simuPOP, numpy, scipy, pprint, subprocess (call), sys, time, glob, fileinput, re, math, os

<b>To Run:</b>
<p>1) download master branch and place all files into a single directory.

2) open a terminal window and cd into the directory containing the scripts.

3) make sure .py and .R scripts are executable (e.g. chmod u+x <file_name.py>).

4) run wrapper scripts locally (e.g. ./simuHybrid_wrapper_disperseArchi.R)</p>

  - there are three wrappers, one for each of three types of genetic architecture: "dispersed", "interspersed", and "modular".
  - running a given wrapper will call the corresponding *.py script and run the simulations. (e.g. "simuHybrid_wrapper_intersperseArchi.R" will call "simuHubrid_schumer_model_intersperseArchi.py")
 
<b>Output:</b>
<p>Running a single .R wrapper will generate directories within the same directory containing the .R and .py scripts (one for each population sized simulated; this is specified by a vector of parameter values in the .R wrappers).
Within each of these directories will be .txt files that contain allele frequencies at each simulated locus, recorded every 100 generations for a total of 1000 generations.</p>

  - The header line of the .txt output files specifies the generation (first column) and the locus for which allele frequencies are recoreded (following columns; one column for each locus).
  - Every 11 rows after the header line contain allele frequencies reported for a single simulated hybrid swarm.
  - There are therefore (n*11)+1 lines in each file, where n is the number of independent hybrid swarms that were simulated with a given combination of paramter values.
  - A seperate .txt file in generated for every different combination of parameter values hybrid swarms are simulated under.
  - The parameter values are specified in the file name. One should check that the parameter values reported in the file names match those specified in the .R wrappers.
