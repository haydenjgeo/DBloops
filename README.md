# DBloops
Matlab functions and scripts for obtaining grain size distribution of coarse (6.3 cm+) particles in nonplanar deposits of granular material

This repository contains all necessary Matlab code to recreate the analyses presented in Jacobson et al., (2023). This includes the clustering code DBloops, the machine learning point cloud descriptor calculation tool Terpukto (Weidner, 2023), the clustering code G3Point (Steer, 2023), and the function confusionmatStats (Cheong, 2023). 

To recreate the analyses presented in Jacobson et al., (2023), users must first download the 'GSD code' folder from this repository. 

Next, users must download the point cloud file 'TrainingRegions.txt' from the USGS data release found at https://doi.org/10.5066/P9QQ0AR2 and place it in the main 'GSD code' folder. This point cloud (3.3 GB) will be copied, subdivided, and modified during recreation of analyses, so ensure a minumum of 14 GB of storage are available before running code. 

Finally, users must open the folder 'ReproduceResults', open the script 'MasterScript.m', and run the script from the main 'GSD code' folder. Note that as precalculated point cloud descriptors (features) are included as the file 'Feat3.mat' it is not necessary to run the script 'SubregionFeatures.m' which takes several hours to run, so the command to execute this script is commented out in 'MasterScript.m'. Users may uncomment and run this code if desired. 

Descriptions of all functions and scripts called by MasterScript.m are included as comments in the script.  

Note on minor edits made to Terpunkto and G3Point: 

The version of Terpunkto uploaded here is renamed 'Terpunkto2scan2geo.m' due to the addition of lidar scan index and geometric roughness based features. Minor edits to G3Point have been made, including commenting out the plotting code to save time - visualization is critical for G3Point optimization, so do not use this version of G3Point independently. Use the version found at the link in the references below. The 'param.csv' file in G3Point has been modified with the best parameters identified from trial and error testing. Default parameters are saved in 'param - backup.csv'. The path name requirements in 'G3Point.m', 'loadptCloud.m', and 'defineparmaeters.m' have also been removed so G3Point can be run from the main 'GSD code' folder. 




References

Cheong, A. (2023). confusionmatStats, Matlab Central File Exhange, https://www.mathworks.com/matlabcentral/fileexchange/46035-confusionmatstats-group-grouphat

Jacobson, H., Walton, G., Rengers., K., and Barnhart, K., (2023). A method to Obtain Remotely Sensed Grain Size Distributions from Nonplanar Granular Deposits. In prep. 

Steer, P. (2023). G3Point, Github Repository,  https://github.com/philippesteer/G3Point

Weidner, L. (2023). Terpunkto, Github Repository, https://github.com/lmweidner/terpunkto
