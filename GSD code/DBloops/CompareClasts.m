function rmse = CompareClasts(cloud1, scsm, sc, trucloud,plt)
% Description: This function downsamples the input 'cloud' to a rough value of 'dsv' points per square meter
%
% Assumptions: cloud1 and trucloud must be formatted exactly as stated below 
%
% Inputs:
%   cloud1 - resulting cloud from DBloops, subsampled [x,y,z,clast index 1, clast index 2, original point index matching trucloud] 
%   scsm - single clasts from decomposed multiclasts (cloud 1 clast index 1)
%   sc - big epsilon single clasts (cloud 1 clast index 2)
%   trucloud - indexed full density clast cloud [X Y Z clast_index]
%   plot - 0 for no plot, 1 for plot 
% 
% Outputs:
%   rmse - based on error in X-direction for each clast 
%
% Revision: R2022a
% Author: Hayden Jacobson
% Date: June 28 2023

%Big to do: When FP points are present, they have no clast index in the
%true cloud. Need to deal with these ASAP. 

gsdscsm = []; %initiate vector of sizes for scsm clasts
gsdscsmtru = []; %vector of "true" sizes to be directly compared to scsm clasts
for k = 1 : length(scsm) % # of scsm clasts
    ind =  find(cloud1(:,4) == scsm(k)); %this is the clast index variable
    tempcloud = cloud1(ind,1:3); %clasts defined by scsm 
    [coeff, score, latent] = pca(tempcloud); %PCA
    trans = tempcloud*coeff; %get size from PCA
    if size(trans,1)<3 %defined by too few points - unable to get size so set to 0 
        gsdscsm(k)=0;
    else
    gsdscsm(k) = max(trans(:,2))-min(trans(:,2)); %define size w/ PCA2
    end
    pointinds = cloud1(ind,6); %original (true) clast indexes for testing cloud clast defined by scsm
    tcind2=mode(trucloud(pointinds,4)); %most common original index, most similar clast 
    ind3=find(trucloud(:,4) == tcind2); %get points with that index from the full density cloud 
    if length(ind3)~=0 %if clast object exists
        tempcloud = trucloud(ind3,1:3); %use PCA to define size again 
        [coeff, score, latent] = pca(tempcloud); 
        trans = tempcloud*coeff;
%         gsdscsmtru(i) = max(trans(:,2))-min(trans(:,2));
          gsdscsmtru(k) = prctile(trans(:,2),99)-prctile(trans(:,2),1);
    else  gsdscsmtru(k)=0; 
    end
end

gsdsc = []; 
gsdsctru = [];
for r = 1 : length(sc) % # of sc clasts
    ind =  find(cloud1(:,5) == sc(r)); 
    tempcloud = cloud1(ind,1:3); 
    [coeff, score, latent] = pca(tempcloud); 
    trans = tempcloud*coeff;
    if size(trans,1)<3
        gsdsc(r)=0;
    else
    gsdsc(r) = max(trans(:,2))-min(trans(:,2));
    end
    tcind = [];
    pointinds = cloud1(ind,6);
    tcind2=mode(trucloud(pointinds,4));
    ind3=find(trucloud(:,4) == tcind2);
    if length(ind3)~=0
        tempcloud = trucloud(ind3,1:3);
        [coeff, score, latent] = pca(tempcloud); 
        trans = tempcloud*coeff;
%         gsdsctru(i) = max(trans(:,2))-min(trans(:,2));
        gsdsctru(r) = prctile(trans(:,2),99)-prctile(trans(:,2),1);
    else gsdsctru(r)=0;
    end
end

%error distribution for clasts: 
gsdest=[gsdsc gsdscsm]; % singles and decomposed singles 
gsdtru=[gsdsctru gsdscsmtru];

%option to identify false clasts composed of outliers and remove them from error calculation
% fcount=0;
% for i=1:length(gsdtru)
%      if gsdtru(i)==0
%          fcount=fcount+1;
%          gsdtru(i)=NaN;
%          gsdest(i)=NaN;
%      end
% end
% sprintf('number of false clasts: %u',fcount)

resids=gsdtru-gsdest; 
rmse = sqrt(mean(resids.^2,'omitnan'));

%plot
if plt==1
figure
%error will be in units of m, bins desired every cm (0.01)
edges = [floor(min(resids)*100)/100:0.01:ceil(max(resids)*100)/100];
histogram(resids,edges);
xlabel('Grain Size Error (m)','fontsize',16)
ylabel('Number of Grains','fontsize',16)
% title(sprintf('GSD comparison: RMSE (m) = %0.5f',rmse))

figure
scatter(gsdest,gsdtru,3,'k','filled')
hold on
dim=max([gsdest gsdtru]);
plot([0 dim+.2], [0 dim+.2])
xlim([0 dim+.2])
ylim([0 dim+.2])
xlabel('Estimated Grain Size (m)','fontsize',16)
ylabel('True Grain Size (m)','fontsize',16)
% title(sprintf('GSD comparison: RMSE (m,x-dir) = %0.5f',rmse))

% figure - trying to bound error percentiles  
% data(:,1)=gsdest;
% data(:,2)=gsdtru;
% plot(data(:,1),data(:,2),'.k','MarkerSize',5)
% hold all
% [N,c] = hist3(data(:,1:2),[20 20]);
% x = c{1};
% y = c{2};
% [X Y] = meshgrid(x,y);
% contour(X,Y,N',[5 10 50 100],'LineWidth',2)
% hold all
% plot(linspace(min(data(:,1)),max(data(:,1))),linspace(min(data(:,1)),max(data(:,1))),'-k','LineWidth',3)
% plot(0.0035*ones(1,100),linspace(min(data(:,2)),max(data(:,2))),'--r','LineWidth',3)
% plot(linspace(min(data(:,1)),max(data(:,1))),0.0035*ones(1,100),'--r','LineWidth',3)


end


end


