%% Script to recreate desired clouds after filtering out clasts identified as undersized 
%Start with full cloud: I already know that field 7 is the clast index and
%9 is the classification index, so these values were manually set below 

indexedcloud=readmatrix('PlaceholderCloud.txt'); %cloud name 
% get GSD of clast points
idcloud=indexedcloud(:,[1 2 3 7]);
[gsd, vals] = GetGSD(idcloud,99.5);
% find undersized clasts 
id=find(gsd<0.063);
length(id)
small2 = vals(id); %small clast indexes from original cloud
tbc=ismember(indexedcloud(:,7),small2);
%change classification of those clast points to 1 (matrix) and set index to NaN
idfield=indexedcloud(:,7);
classfield=indexedcloud(:,9);
idfield(tbc==1)=NaN; %matrix index is NaN
classfield(tbc==1)=1; %matrix classification is 1
indexedcloud(:,7)=idfield;
indexedcloud(:,9)=classfield;
writematrix(indexedcloud,'FH_NoUndersized_AllFields.txt');

%% verify removal 
indexedcloud=readmatrix('FH_NoUndersized_AllFields.txt'); 
idcloud=indexedcloud(:,[1 2 3 7]);
[gsd, vals] = GetGSD(idcloud,99.5);
id=find(gsd<0.063); 
length(id)
%% compare distributions from use of different percentiles 
gsdtrue1 = GetGSD(indexedcloud,99); %get GSD from indexed clast cloud 
[idx2, gsdtrue2] = GetGSD(indexedcloud,100);
figure
[f1,x1]=ecdf(gsdtrue1); %y is size, x is percentile
[f2,x2]=ecdf(gsdtrue2);
plot(x1,f1,'LineWidth',2)
hold on
plot(x2,f2,'LineWidth',2)
legend({'1-99%','0-100%'},'fontsize',18)


%% looking at individual clasts if desired
%Have to drop one index because clast indexes started at 0 rather than 1
cin=15;
cin=cin-1;
newclast=indexedcloud(indexedcloud(:,4)==cin,1:3);
in=newclast;
[coeff, score, latent] = pca(in(:,1:3));
trans = in*coeff;
gs = max(trans(:,2))-min(trans(:,2));
figure %original orientation
scatter3(newclast(:,1),newclast(:,2),newclast(:,3))
figure %translated
scatter3(trans(:,1),trans(:,2),trans(:,3))
