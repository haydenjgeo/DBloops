% Description: This script runs the necessary functions to prepare data for
% DBloops and then applies DBloops to the manually delineated clast points
% only 
% 
% Outputs: DBloops figure displaying estimated and true GSD
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%%
% fcloud=readmatrix('FH_NoUndersized_AllFields.txt');%load in full point cloud 
fcloud=readmatrix('TrainingRegion.txt');
ccloud=fcloud(~isnan(fcloud(:,7)),[1 2 3 7]);%select clast points with clast index, non-clast points have NaN index
c2cloud=ccloud(:,[1 2 3]);%cloud of clast points no index
gsdtrue = GetGSD(ccloud,99.5); %get GSD from indexed clast cloud 
dens=1500; %downsampling density 
[clouda, idx] = DownCloud(c2cloud,dens); %downsampled cloud
np=18; %number of points (corresponds to 'dens' of 1500 ppm^2)
esfa=2.6; esfb=0.6; %epsilon scaling factors; optimization performed w/ these values
%%
[ea, eb] = EstimateEpsilon(np, esfa, esfb, clouda); %get scaled epsilon values
[gsd, cloudb, scsm, sc, prmse, frmse]=DBloops(gsdtrue,clouda,np,ea,eb,1); %cluster with DBloops
gsdlength=length(gsd);
%% save error estimates and number of clasts, print figure 
save F9_ClusterTest.mat prmse frmse gsdlength
print('F9_ClusterTest', '-dpng', '-r400') %save figure
