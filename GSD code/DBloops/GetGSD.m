function [gsd, vals] = GetGSD(indexedcloud,pr)
% Description: This function obtains the grain size distribution of a point cloud of indexed clasts
% 
% Assumptions: The cloud only contains only [X Y Z clast_index] 
%   If any matrix present, index must be NaN for correct removal 
%
% Inputs:
%   indexedcloud - input cloud formatted as matrix 
%   pr - percentile range of points aligned to second eigenvector for size estimation; 
%       99 indicates 1-99, 100 indicates 0-100
% 
% Outputs:
%   gsd - grain size distribution of clasts 
%
% Revision: R2022a
% Author: Hayden Jacobson
% Date: June 28 2023



% gsd = zeros(1,length(unique(indexedcloud(:,4)))); %preallocate GSD vector
vals = unique(indexedcloud(:,4)); %identify all clast indexes 
vals=vals(~isnan(vals)); %if cloud of matrix and clasts is passed remove nans 
gsd=zeros(1,length(vals));

for i = 1 : length(gsd)
    valind=find(indexedcloud(:,4)==vals(i));
    in = indexedcloud(valind,1:3);
    [coeff, score, latent] = pca(in);
    trans = in*coeff;
    if pr==99
        gsd(i) = prctile(trans(:,2),99)-prctile(trans(:,2),1);
    elseif pr==99.5
        gsd(i) = prctile(trans(:,2),99.5)-prctile(trans(:,2),0.5);
    elseif pr==100
        gsd(i) = max(trans(:,2))-min(trans(:,2));
    end
end
end
