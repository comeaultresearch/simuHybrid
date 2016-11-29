#!/usr/bin/env Rscript
#./simuHybrid_wrapper_intersperseArchi.R

#################################################################
# wrapper to run "simuHybrid_schumer_model_intersperseArchi.py" #
# Aaron A. Comeault (aacomeault@gmail.com)                      #
# 18 Nov. 2016                                                  #
#################################################################
rm(list=ls())

# specify parameters to use during simulations:
ns = c(50, 100, 1000)
s_dmis = c(0.000, 0.001, 0.01, 0.05, 0.1)
s_adds = seq(0.0, 0.1, 0.01)
h = .5
rs = c(0.5, 0.1, 0.05) # .5 is 40 cM region, 0.1 is 8 cM region, 0.05 is 4 cM region, 0.01 is 0.8 cM.


# call the "simuHybrid_schumer_model_intersperseArchi.py" script
# for each parameter combination:
for (n in ns) {
  for (rec.rate in rs) {
    for (s_dmi in s_dmis) {
      for (s_add in s_adds) {
        
        # call the .py script to run simulations:
        system(paste("./simuHybrid_schumer_model_intersperseArchi.py",  s_dmi, s_add, h, rec.rate, n, sep=" "))
        
      }
    }
  }
}
