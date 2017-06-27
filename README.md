# LRd-code
Matlab code associated with LRd model, for REU 2017 project

fixed-point-revision branch (LMM) 

fixedpt_search: New script, loops through bcls and tries to identify a fixed point for each one. Fixed points are stored in a 17-times-number-of-bcls matrix named allfp. 

lrd_p2p: Made changes to make it compatible with nsoli.m, in addition to pacedown.m. It should still be compatible with pacedown.m, but I haven't tested this. The problem is that nsoli expects that the function that we're evaluating to only have one input (the initial condition), so we can either modify lrd_p2p or nsoli, so I chose the former. Modified lrd_p2p to read in any missing arguments after yin (bcl, ncyc, subdiv_per_cyc) from a file, if the file is available, before reverting to default values. Changed 'lrdparams' filename to 'lrdinputs', since the former was used for different purposes than are needed here. Folder for saving data is now defined at the beginning.

lrddata_1cell_b1000.mat: The final state of Y in this file is the fixed point that I found for bcl = 1000 on my Surface Pro 3. In order to get the fixed point to be loaded properly by fixedpt_search.m, the lrddata...mat file must be moved into the lrddata/ directory. 
