function [gsd, cloud1, scsm, sc, partial_rmse, full_rmse] = DBloops(gsdtrue, cloud, np, ea, eb, plt)
%
% Description: This function clusters point cloud data of granular deposits
% to obtain a GSD of particles >6.3 cm in size
%
% Assumptions: Imput data consists of only clast points, XYZ are first
%
% Inputs:
%   gsdtrue - vector of true GSD, likely obtained from GetGSD.m (can also
%   use false data e.g. [0 0 0] if not interested in rmse or true GSD not
%   available)
%   cloud - input clast cloud formatted as matrix [X Y Z (other SFS)]
%   np - number of points 
%   ea - small epsilon value (lots of outliers) 
%   eb - big epsilon value (lots of multiclasts)
%   plt - 0 for no plot, 1 for standard plot, 2 for simplified tiledlayout plot
% 
% Outputs: Variables 
%   gsd = gsd returned by DBloops (m)
%   cloud1 - [X Y Z scsm sc] returns XYZ with identified clast indexes
%       Note: outliers have index of -1
%   scsm - single clast indexes from decomposed multiclasts (cloud 1 clast index 1)
%   sc - big epsilon single clast indexes (cloud 1 clast index 2)
%   partial_rmse - rmse fropm 5-95% of distribution (m)
%   full_rmse - rmse from full distributions (m)
% Outputs: Command Window
%   outliersa - # ea outliers before outlier assignment 
%   outliersb - # eb outliers before outlier assignment 
%   outliersa2 - # ea outliers after outlier assignment 
%   outliersb2 - # eb outliers after outlier assignment 
%   KS2 test results
%
%  
% Note: variables with "ea" or "a" generally refer to small epsilon (lots of outliers),
% similar for "eb" or "b" referring to large epsilon (lots of multiclasts)
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

cloud = cloud(:,1:3); % load cloud
cloud1 = cloud;

% Run DBSCAN 2x 
[idxa] = dbscan(cloud, ea, np); 
cloud1 = [cloud1 idxa]; %index 2 col 4
[idxb] = dbscan(cloud, eb, np); 
cloud1 = [cloud1 idxb]; %index 3 col 5

%check number of identified clasts fpr each epsilon
numa = length(unique(idxa))-1; %index -1 represents outliers, drop one index
numb = length(unique(idxb))-1;

%find and print # outliers for each epsilon
outliersa = length(idxa(idxa==-1))
outliersb = length(idxa(idxb==-1))

% Check if eb clasts are multiclasts by comparing sizes
%This relies on the low epsilon value not oversplitting clasts
mc = []; %multiclast eb clast indexes (note this is not a continuous variable)
sc = []; %single clast eb clast indexes (note this is not a continuous variable)
outr = 1+(outliersb/outliersa); %outlier ratio used to define multiclasts
for i = 1:numa %number of unique ea (small epsilon) clasts
    tempc=cloud1(idxa==i,:); %points for relevant ea clast
    t = mode(tempc(:,5)); %index for eb clast most similar to ea clast
    if (length(idxb(idxb==t))) < outr*(length(idxa(idxa==i))) %similarity of size (Equation 1 in main text)
        sc = [sc t];%eb cluster probably single clast, add to single clast list
    else
        mc = [mc t];%add eb cluster index to probable multiclast list
    end 
end
sc=unique(sc); %collapse indexes in sc and mc vectors 
mc=unique(mc);
sc = setdiff(sc, mc); %return data in sc that is NOT in mc due to potential overlap
scx = setdiff(unique(idxb),[sc mc]); %clasts in idxb but not sc/mc - generally more outlying clasts that were assigned outliers at ea
scx(scx==-1)=[]; %remove outliers
sc = [sc scx']; %combine lists 

%Outlier assignment step 1: assign all large eb outliers to nearest eb clast if within 10 cm 
%for single clasts only
outs = find(idxb==-1); %outlier indexes
cloud2 = cloud1(find(idxb~=-1),:); %make cloud of non-outliers 
for i=1:length(outs) %look at every outlier (hi ep)
    x = cloud(outs(i),:); %outlier observations   
    d = pdist2(x,cloud2(:,1:3)); %outlier vs cloud of non outliers point-to-point distances
    [m,r] = min(d); %r is index of closest non-outlier to point
    if m<0.1 %within 10 cm 
        cloud1(outs(i),5)=cloud2(r,5); %give outlier point correct eb clast index
    end
end
outliersb2 = length(find(cloud1(:,5)==-1)) %show reduced number of eb outliers

% Outlier assignmnet step 2: assign all ea outliers to nearest ea clast (for points in eb multiclasts only)
for i = 1:length(mc)
    ins = cloud1((find(cloud1(:,5)==mc(i) & cloud1(:,4)~=-1)),:);%indexes of eb multiclast points that are NOT ea outliers
    outs1 = find(cloud1(:,5)==mc(i) & cloud1(:,4)==-1); %indexes of eb multiclast points that ARE ea outliers
    for j=1:length(outs1) %for each outlier...
        x=cloud(outs1(j),:); %make a one point cloud for the outlier...
        d=pdist2(x,ins(:,1:3)); %and compare outlier point to rest of points in MC
        [m,r]=min(d); %find closest distance and index
        if m<0.1 %if within 10 cm 
            cloud1(outs1(j),4)=ins(r,4); %assign index of nearest lo ep cluster to outlier
        end
    end
end
outliersa2 = length(find(cloud1(:,4)==-1)) %show reduced number of ea outliers

% Create list of all ea single clasts that are part of eb multiclasts
scsm = []; %eb single clasts 
for i = 1:length(mc)
    r = find(cloud1(:,5)==mc(i)); %eb clasts that are on multiclast list 
    tlist = unique(cloud1(r,4)); %ea clast indexes that compose relevant multiclast
    scsm = [scsm tlist']; %ea clast list
end
scsm = unique(scsm); %collapse ea clast list
scsm(scsm==-1)=[]; %remove outlier index

% Get GSD of ea clasts that are part of eb multiclasts 
gsdscsm = []; % storage for GSD
smallscsm = 0; % counter for clasts with <6 points 
for i = 1 : length(scsm) % #clasts
    ind =  find(cloud1(:,4) == scsm(i)); %relevant clast index
    tempcloud = cloud(ind,:); %temporary point cloud 
    [coeff] = pca(tempcloud); 
    trans = tempcloud*coeff; 
    if length(trans(1,:))>=2
        gsdscsm(i) = max(trans(:,2))-min(trans(:,2));
    else
        smallscsm=smallscsm+1;
    end
end

%Get GSD of eb single clasts
gsdsc = []; 
smallsc = 0; % counter for clasts with <6 points
for i = 1 : length(sc) %#clasts
    ind =  find(cloud1(:,5) == sc(i)); %relevant clast index
    tempcloud = cloud(ind,:); %temporary point cloud
    [coeff] = pca(tempcloud); 
    trans = tempcloud*coeff;
    if length(trans(1,:))>=2
        gsdsc(i) = max(trans(:,2))-min(trans(:,2));
    else
        smallsc=smallsc+1;
        %clusters with <6 points have insufficient dimensnlty
    end
end

%Merge eb and ea single clast GSDs to get final GSD
gsd=[gsdsc gsdscsm];
if isempty(gsd)==0 % if gsd is not empty 
    [f6,x6]=ecdf([gsdsc gsdscsm]);% create cdf
end
if x6(1)==0 && f6(1)==0
    x6(1)=[]; f6(1)=[]; %remove false 0 that may arise from ecdf
end

%Make CDF of true GSD for plotting/error calculation
[f1,x1] = ecdf(gsdtrue); %ecdf for true data


%use linspace to calculate residuals and RMSE
[xf1 I] = unique(x1); %true
ff1 = f1(I);
[xf6 I] = unique(x6); %estimated 
ff6 = f6(I);

%blow up error if <2 clasts recognized and don't try to calculate residuals or plot
if length(xf6)<2 
    rmse = 100;
else
intax=interp1(ff1,xf1,0:.01:1);
intbx=interp1(ff6,xf6,0:.01:1);

%calculate residuals
residsx=intbx-intax;
full_rmse = sqrt(mean(residsx.^2,'omitnan')) %print full rmse

%filter residuals to 5th to 95th percentile of distribution 
intax2=interp1(ff1,xf1,.05:.01:0.95);
intbx2=interp1(ff6, xf6,.05:.01:0.95);
residsx2=intbx2-intax2;
partial_rmse = sqrt(mean(residsx2.^2,'omitnan')) %print partial rmse

%plot and format 
if plt==1 %standard figure 
    ff=figure;
    ff.Position=[0 0 900 400];
    stairs(x1,f1,'color','k','LineWidth',3);
    hold on
    stairs(x6,f6,'LineWidth',3,'color','b')
    plot(abs(residsx),0:.01:1,'color',[0.4660 0.6740 0.1880],'LineWidth',2);
    settings.matlab.fonts = 'Tahoma';
    if max(gsdtrue)>max([gsdsc gsdscsm])
        set(gca,'yscale','linear','xscale','linear','xlim',[0 max(gsdtrue)],'fontname','helvetica');
    elseif max(gsdtrue)<max([gsdsc gsdscsm])
        set(gca,'yscale','linear','xscale','linear','xlim',[0 max([gsdsc gsdscsm])],'fontname','helvetica');
    end
    xlabel('Grain Diameter (m)','fontsize',16,'fontname','helvetica');
    ylabel('Proportion Passing','fontsize',16,'fontname','helvetica');
    legend('Training Data','DBloops Result','|Residuals|','fontsize',14,'location','southeast','fontname','helvetica');
end

if plt==2 %tilledlayout (used in F11_CrossValidation.m) 
    nexttile
    stairs(x1,f1,'color','k','LineWidth',2);
    hold on
    stairs(x6,f6,'LineWidth',2,'color','b')
    plot(abs(residsx),0:.01:1,'color',[0.4660 0.6740 0.1880],'LineWidth',2);
    if max(gsdtrue)>max([gsdsc gsdscsm])
        set(gca,'yscale','linear','xscale','linear','xlim',[0 max(gsdtrue)]);
    elseif max(gsdtrue)<max([gsdsc gsdscsm])
        set(gca,'yscale','linear','xscale','linear','xlim',[0 max([gsdsc gsdscsm])]);
    end
end

[h3,p3]=kstest2(gsdtrue, [gsdsc gsdscsm]);
formspec = 'KS2 test: h=%u, p=%f';
ks2string = sprintf(formspec,h3,p3) %print KS2 test results 
end
end