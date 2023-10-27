function [cloud] = InterpLabel(FullCloud,LabeledCloud)
% Description: This function interpolates one scalar field from a subsampled cloud to a  full density cloud using KNN (N=2)
%
% Assumptions: Clouds overlap, first three indexes of input clouds are X,Y,Z, and labeled cloud only has 1 SF
%
% Inputs:
%   FullCloud - Full density cloud with any # of scalar fields 
%   LabeledCloud - Cloud with 1 scalar field to be interpolated [X Y Z label]
% 
% Outputs:
%   cloud - full density cloud with interpolated scalar field 
%
% Revision: R2022a
% Author: Hayden Jacobson
% Date: June 28 2023

FullCloud = FullCloud(:,1:3);
Labtemp=LabeledCloud(:,1:3); %just use first three rows for KNN

idx = knnsearch(Labtemp,FullCloud,'K',1); %get indexes of nearest labeled points for each point in full cloud
 
for i = 1:length(FullCloud)
    FullCloud(i,4)=LabeledCloud(idx(i),4); %apply label of nearest point
end

cloud=FullCloud;

end

