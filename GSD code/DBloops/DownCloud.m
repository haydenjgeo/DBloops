function [cloud1, idx, dens] = DownCloud(cloud, dsv)
% Description: This function downsamples the input 'cloud' to a rough value of 'dsv' points per square meter
%
% Assumptions: First three indexes of input cloud are X,Y,Z
%
% Inputs:
%   cloud - input cloud formatted as matrix 
%   dsv - desired density    
% 
% Outputs:
%   cloud1 - downsampled cloud
%   idx - index of subsampled points in original cloud
%   desn - density of cloud
%
% Revision: R2022a
% Author: Hayden Jacobson
% Date: June 28 2023


%Quick PCA analysis to estimate surface area of cloud  
[coeff,score,latent] = pca(cloud(:,1:3));
trans = cloud(:,1:3)*coeff;
PC1= prctile(trans(:,1),99)-prctile(trans(:,1),1);
PC2= prctile(trans(:,2),99)-prctile(trans(:,2),1);
area = PC2*PC1*0.85; 
%scale area by 85% to account for deposit curvature - based on experience

%full density estimate 
dens=length(cloud)/area;

%find appropriate number of points for downsampling 
dspts=floor(area*dsv);

% downsample cloud and return sample indexes 
[cloud1, idx] = datasample(cloud,dspts);
idx=idx';

end
