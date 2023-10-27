% Description: This script conducts manual (trial and error) optimization
% of DBloops parameters by testing DBloops over a range of MinPts values at different densities
% to identify a minimum density threshold at which both error and epsilon begin to
% stabilize and determine appropriate MinPts values
%
% Outputs: Plots showing relationships between density and error, epsilon,
% and MinPts; data for plotting also saved in A567_NumPointsOptimization.mat
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%% Load point clouds
%Read in indexed and unindexed clouds containing only clasts (indexed is "true")
numareas=5;
for i=1:(numareas)
    fcloud=readmatrix(sprintf('ab%d',i)); %full cloud
%clast points with index; I know that col 7 is clast index, 9 is classification
    ccloud=fcloud(~isnan(fcloud(:,7)),[1 2 3 7]);
    tc{i,1}=ccloud;
    [gsdtrue, vals] = GetGSD(ccloud,99.5);
    gsdt{i,1}=gsdtrue;
    numrock(i)=length(gsdt{i});
    c2cloud=ccloud(:,[1 2 3]);
    c1{i,1}=c2cloud;
end

%% Select density values to be tested
density=[2500 2250 2000 1750 1500 1250 1000 900 800 700 600 500 400 300 200 100 50]; %no downsampling for first test if doing full density..so slow dont do it lol

%% Run brute force optimization (minimize CompareClasts error) 
for q=1:numareas %#full density clouds (1:5)
for j=1:length(density) %different cloud densities (1:17)
    dens=density(j);
    [clouda, idx] = DownCloud(c1{q,1},dens);
    nps=[3:26]; %estimated range of appropriate values

    %fill all fields with NaNs at numpts less than tested (here just 1-2)
for w=1:(min(nps)-1)
sc2(w)=NaN;
scsm2(w)=NaN;
eps(w)=NaN;
error(w)=NaN;
rmse2(w)=NaN;
end

for k=min(nps):max(nps)
    np=k;
    %Note DBscan will error out for unreasonable point guesses, look for length [scsm+sm]=0 
    esfa=2.6; esfb=0.6; %Test results only valid for this scaling ratio
    [ea, eb] = EstimateEpsilon(np, esfa, esfb, clouda);
    [gsd, cloudb, scsm, sc, rmse(k), fullrmse]=DBloops(gsdt{q},clouda,np,ea,eb,0);
    scsm2(k)=length(scsm);
    sc2(k)=length(sc);
    cloudb = [cloudb idx];
    clear i;
    if length([scsm sc])>2
        rmse2(k)=CompareClasts(cloudb, scsm ,sc, tc{q,1},0);
    else %blow up error if no clasts recognized
        rmse2(k)=100000;
    end
    eps(k)=ea*esfa; %EstimateEpsilon does not provide base epsilon
end

error=sqrt(rmse2.*rmse); %geometric mean - error with units of meters

[val2,idx2]=min(error);
erroropt(q,j)=val2;
epsilonopt(q,j)=eps(idx2);
npopt(q,j)=idx2; %after adding NaN population, idx2 is the numpoints at min error 

end
end

save A567_NumPointsOptimization.mat density epsilonopt npopt erroropt numareas

%% load data if desired
load('A567_NumPointsOptimization.mat')

%% Figure A5: optimized numpoints at different densities
numareas=5;
ff=figure;
ff.Position=[0 0 700 500];
colors=({'r' 'g' 'b' 'm' 'k'});
for t=1:numareas %numsites; one curve for each 
    scatter(density,npopt(t,:),50,char(colors(t)),'filled')%,'color',([1 1 1]./t))
    hold on
end
ddensity=repmat(density,5,1);
f=fit(ddensity(:),npopt(:),'poly2');
p=plot(f);
p.LineWidth=2;
p.Color='black';
xlabel('Clast Point Density (points/m^2)')
ylabel('Optimized MinPts')
lgd=legend({'1' '2' "3" "4" "5",'2^{nd} Order Polynomial'},'location','southeast');
lgd.Title.String='Training Data';
print('A5_MinPointsVDensity', '-dpng', '-r400')


%% Figure A6: Clast size error at different densities
numareas=5;
ff=figure;
ff.Position=[0 0 700 500];
colors=({'r' 'g' 'b' 'm' 'k'});
for t=1:numareas
    hold on
    scatter(density,erroropt(t,:),50,char(colors(t)),'filled');

end
xlabel('Clast Point Density (points/m^2)')
ylabel('Clustering Error (m)')
lgd=legend({'1' '2' "3" "4" "5"},'location','northeast');
lgd.Title.String='Training Data';
print('A6_ErrorVDensity', '-dpng', '-r400')

%% Figure A7: optimized epsilon w different densities 
ff=figure;
ff.Position=[0 0 700 500];
numareas=5;
colors=({'r' 'g' 'b' 'm' 'k'});
for t=1:numareas %numsites; one curve for each 
    scatter(density,epsilonopt(t,:),50,char(colors(t)),'filled')%deal with ,'color',([1 1 1]./t))
    hold on
end
xlabel('Clast Point Density (points/m^2)')
ylabel('Optimized Epsilon (cm)')
lgd=legend({'1' '2' "3" "4" "5"},'location','northeast');
lgd.Title.String='Training Data';
print('A7_EpsilonVDensity', '-dpng', '-r400')
