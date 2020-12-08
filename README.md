# NeuroFdpToolbox
MATLAB toolbox to automatically detect and analyse fractional Lévy motions of propagating neural activity patterns, developed by Dr. Pulin Gong's group at University of Sydney. This toolbox includes the generation of a spiking neural network model, and the detailed processes for analysing our simulated and experimental data.
If you use our code in your research, please cite us as follows:

Liu Y., Long X., Martin PR, Solomon SG and Gong P., Lévy walk dynamics explain gamma burst patterns in primate cerebral cortex. Under Review, 2020. 

## Network generation
A spiking neural circuit simulation model edited by Yifan Gu, Yuxi Liu, James Henderson, Guozhang Chen. The instruction is detailed in the readme file in the sub-folder.


## Simulation data analysis
The matlab function for generating the neuron microcircuit is main_Gamma.m 
The detection, analysis and visualization of the fractional propagating patterns in the simulation data can be found in the folder model_data_analysis. 
The code to generate the figure 4-6 in the paper "Lévy walk dynamics explain gamma burst patterns in primate cerebral cortex" can be found in the sub-folder experimental_data_analysis/New/GammaPaperFig1-4.

## Experimental data analysis
The detection, analysis and visualization of the fractional propagating patterns in the experimental data can be found in the main function Project1.m in the sub-folder experimental_data_analysis/Toolbox_CSC. An example movie of the experimental data can be found:

![Example superdiffusive gamma burst pattern movie](https://github.com/longxian319/PhD_XL/blob/master/GammaDynaPatt/example%20movies/GammaBurstPatterns.avi)

The full data can be shared upon requested. The code to generate the figure 2 and figure 3 in the paper "Lévy walk dynamics explain gamma burst patterns in primate cerebral cortex" can be found in the sub-folder experimental_data_analysis/Toolbox_CSC/p1.


## Authors

* **Yuxi Liu** - *Model generation and analysis* - [yliu2521](https://github.com/yliu2521)
* **Xian Long** - *Experimental data analysis* - [Xian Long](https://github.com/longxian319)
* **Pulin Gong** - *Coodinator* - pulin.gong@sydney.edu.au
