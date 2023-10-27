% Description: This script uses preexisting point cloud indexes to recreate
% the manual division of the training region into 5 subregions (a#.txt). This script
% also checks for and removes any undersized clasts in the training data,
% creating both subregions (ab#.txt files) and the full training region
% (FH_NoUndersized_AllFields.txt) with no undersized clasts. These files
% are used in all analyses. 
%
% Outputs: point clouds of the subregions with undersized clasts included (5) 
% and the full training region (1) and the subregions(5)with no
% undersized clasts
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%% Load full point cloud and subregion point indexes 
%Note cloud contains several slightly undersized (e.g. 6-6.2 cm) clasts when
%GetGSD(cloud,99.5) applied - this script checks for and removes these for later analyses
cloud=readmatrix('TrainingRegion.txt');
load('SubregionIDs.mat');

%% create subregions and then remove any undersized clasts from each
for i=1:5
    indexedcloud=cloud(id{i},:); %indexes for all points in manually divided subregions
    writematrix(indexedcloud,sprintf('Data\\Clouds\\Subregions\\a%d.txt',i));
% get GSD of clast points in subregion
    idcloud=indexedcloud(:,[1 2 3 7]);
    [gsd, vals] = GetGSD(idcloud,99.5);
% find undersized clasts if any exist in original subregions clouds (a1-a5)
    idd=find(gsd<0.063);
    small2 = vals(idd); %small clast indexes from original cloud
    tbc=ismember(indexedcloud(:,7),small2);
%change classification of those clast points to 1 (matrix) and set index to NaN
    idfield=indexedcloud(:,7);
    classfield=indexedcloud(:,9);
    idfield(tbc==1)=NaN; %matrix index is NaN
    classfield(tbc==1)=1; %matrix classification is 1
    indexedcloud(:,7)=idfield;
    indexedcloud(:,9)=classfield;
%create subregions point clouds that contain no undersized material 
    writematrix(indexedcloud,sprintf('Data\\Clouds\\Subregions\\ab%d.txt',i));
end

%% verify no undersized clasts are present in modified point clouds (ab 1-5)
for i=1:5
    indexedcloud=readmatrix(sprintf('ab%d.txt',i)); 
    idcloud=indexedcloud(:,[1 2 3 7]);
    [gsd, vals] = GetGSD(idcloud,99.5);
    idd=find(gsd<0.063);
    length(idd)
end

%% merge clouds back together, but with all clasts removed
for i=1:5
    cloud2{i}=readmatrix(sprintf('ab%d.txt',i));
end

FH = [cloud2{1}; cloud2{2}; cloud2{3}; cloud2{4}; cloud2{5}];
writematrix(FH,'Data\\Clouds\\FH_NoUndersized_AllFields.txt');


