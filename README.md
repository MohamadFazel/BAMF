 Bayesian multiple emitter fitting (BAMF) is a software package developed based on the framework
Described in "Bayesian Multiple Emitter Fitting using Reversible Jump Markov Chain Monte Carlo" 
by Mohamadreza Fazel, Michael J. Wester, Hanieh Mazloom-Farsibaf-Farsibaf, Marjolein B.M. Meddens,
Alexander S. Eklund, Thomas Schlichthaerle, Florian Schulender, Ralf Jungmann and Keith A. Lidke; 
Scientific Reports 9 13791 (2019). 

Subdirectories include

@RJ		      MATLAB codes for the class RJ
mex		      mex-functions for RJMCMC
ptx		      single emitter CUDA-code used in estimating the intensity priors
source		      all the C++ and CUDA source codes
Examples	two examples: example scripts for synthetic and experimental data

Running RJMCMC requires a GPU card, CUDA and MATLAB 2015a or higher.  The user
needs to first compile the C++ and CUDA source codes, putting the results into
mex and ptx, by running the Make script in the top-level RJMCMC directory.
This is necessary as the ptx-files have to be compatible with the users version
of MATLAB and the operating system.  NOTE that the CUDA path on the user's
machine must be set correctly at the top of Make.m .

The BAMF method is capable of localizing emitters in very dense and noisy regions 
of SMLM raw data. The code can be run using either the gun or example scripts. To
Open the gui, cd to the BAMF directory and type RJ in the command window. A guy
will then popup containing multiple sections:

1) RJstruct:     containing parameters required by the algorithm.
2) SMAStruct:    containing parameters used by the single emitter code used to find 
                 empirical prior on emitter intensities.
3) P_Burnin:     probabilities of proposing each jump (step in Gibbs) for burnin chain
4) P_Tiral:      probabilities of proposing each jump for post-burnin part of the chain 
5) Parameters:   experimental parameters
6) Threshold:    thresholds used in filtering out coordinates
7) FrameConnect: parameters used in the frame-connection post-processing step
Use the button "Choose Files" to load data and then hit the "run RJMCMC" button to run 
The code.

The code can also be run using scripts. The Examples directory contains two example of 
data that can be analyzed by RJMCMC. The user should cd in the Example directory in his 
or her MATLAB and run the examples from there.  The synthetic data example runs quickly, 
while the experimental data example requires some time to run.
Please Note that the Actin data file was too  large and could not been uploaded here. This data set is included in the software package in nature website. Please download the software package using this link: https://www.nature.com/articles/s41598-019-50232-x#Sec21
