% Description: This script reproduces all data, figures, and appendix
% figures from Jacobson et al. (2023)
%
% Note: To run this code, you must download 'TrainingRegions.txt' from the 
% USGS data release found at https://doi.org/10.5066/P9QQ0AR2 and place it 
% in the main "GSD code" folder
%
% Outputs: All point clouds and features sets required to produce figures, 
% figures will be saved.  
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%% Create subregions and remove any undersized clasts (several minutes)
CreateSubregions

%% Obtain subregion features (4-12 hours to run) 
SubregionFeatures

%% Conduct Recursive Feature Elimination, create figure A4 (<1 hour)
RFE_A4
close all
%% Produce Figures 8-12 from main text and Table 1 results (~1 hour)
% Figure 8 - comparison of training data and Wolman Count
F8_TrainingWolman
close all
% Figure 9 - test of DBloops on clast points from training data
F9_ClusterTest
close all
% Figure 10 - full workflow (random forest classification + DBloops) 
F10_FullWorkflow
close all
% Figure 11 - block cross validation test of full workflow on subregions
F11_CrossValidation
close all
% Figure 12 - compare G3Point and DBloops 
F12_CompareG3DBloops

%% Produce Appendix Figures 1-3 and 5-9 (several hours due to trial and error optimization)
% Figure A1 - Impact of testing/training split on F-Score and Precision
A1_RatioTest
close all
% Figure A2 - Performance of random forest at variable numbers of trees
A2_NumTrees
close all
% Figure A3 - Performance of random forest at variable number of leaves
A3_LeafDepth
close all
% Figure A5, A6, A7 - MinPts, Error, and Epsilon at different DBloops clast point densities 
A567_NumPointsOptimization
close all
% Figure A8 - Impact of decision threshold (or split point) on Precision
A8_SplitPoint
close all
% Figure A9 - Impact of pre-interpolation density on Precision
%Note "interpolation" is replaced with term "extrapolation" in all text
A9_DensityB4Interpolation
close all