function [ea, eb] = estimateEpsilon(np,esfa,esfb,cloud)
% Description: This function estimates epsilon and provides two scaled values for use with DBloops f
%
% Assumptions: Cloud consists of only clasts or objects to be clustered
%
% Inputs:
%   np - number of points for DBscan, DBloops
%   esfa - ea=epsilon/esfa; generally esfa=2.6
%   esfb - eb=epsilon/esfb; generally esfb=0.6
%   cloud - input cloud formatted as matrix 
%   dsv - desired density    
% 
% Outputs:
%   ea - small epsilon value for DBloops
%   eb - large epsilon value for DBloops
%
% Revision: R2022a
% Author: Hayden Jacobson
% Date: June 28 2023

cloud=cloud(:,1:3); %input cloud must be n*3 matrix (X,Y,Z only)

%knn for point cloud for nearest np points 
[nn, d] = knnsearch(cloud,cloud,'K',np); %n is numpts

%sort knn distances low to high 
dm = sort(max(d,[],2))';
x = 1:length(dm);

%downsample-doesn't need to result in integers. >50 pts convexity uncertain
% dm2=dm(1:length(dm)/50:length(dm)); causes int type error 
% x2=x(1:length(dm)/50:length(dm));
dm2=dm(1:int32(length(dm)/50):length(dm));
x2=x(1:int32(length(dm)/50):length(dm));

%determine convexity 
conc = gradient(gradient(dm2));

%find first + value, over 50 pts threshold of 1e-5 estimated. 
%conc>0 doens't work - variability/small inflections 
wind = find(conc>1e-5,1);
inflection = x2(wind);
dm3=dm(inflection:end);

%identify concave up section only 
liny3 = linspace(dm(inflection), max(dm), length(dm3))';
x3 = linspace(x(inflection), max(x), length(dm3))';

%find "knee", or epsilon, at max distance 
[idx , dist] = knnsearch([x3 liny3],[x3 dm3'],'K',1);
[md mindx]= max(dist);
epsilon = dm3(mindx);

%return scaled epsilon values 
ea = epsilon/esfa;
eb = epsilon/esfb;

% below code plots, rarely needed, just uncomment 
% figure; 
% plot(x,dm); 
% hold on
% scatter(x2,dm2); %query points for convexity 
% plot(x3,liny3,'linewidth',3);
% plot(x3,dm3');
% scatter(x3(mindx),epsilon,100,'filled');
% legend()
% formatting
% settings.matlab.fonts = 'Tahoma';
% formatSpeca = 'KNN distances (m), n=%u points';
% pointstring = sprintf(formatSpeca,n);
% ylabel(pointstring,'fontsize',16);
% xlabel('Points','fontsize',16);
% title('Epsilon Identification','fontsize',18);
% formatSpec = 'Epsilon =%0.5fm';
% epsilonstring = sprintf(formatSpec,epsilon);
% legend('Full KNN','Subsampled KNN','+Concavity', 'KNN +Concavity ',epsilonstring,'fontsize',14);

end