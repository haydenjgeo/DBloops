% Description: This script uses a modified version of Terpunkto (Weidner et al.,
% 2019, 2021) to develop features sets for all subregions. This script
% takes several hours to run. 
%
% Outputs: Full feature sets for all subregions, each with 100,000 core
% points (50,000 for both classes of clast and matrix) used in later
% analyses
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023


%% obtain features - this takes several hours to run
for i=1:5
    feats{i}=terpunkto2scan2geo(sprintf('ab%d.txt',i),[1 2 3 4 5 6 11 12 13 9 8 10],[0.3 0.2 0.15 0.1 0.075 0.05 0.03 0.02 0.01],0,1,1,50000,1,0.01,0.3); 
end

save('Terpunkto\\Features\\Feat3.mat','feats')