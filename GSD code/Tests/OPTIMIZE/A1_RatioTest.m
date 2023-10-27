% Description: The following code evaluates the impact of the
% testing/training ratio on both F-Score (clast class) and Precision, where
% F-Score is evaluated at a confidence of 0.5 and Precision at a confidence
% of 0.65.
%
% Outputs: Plots showing relationship between training/testing proportion,
% F-Score, and Precision
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%% load precalculated features
load('Feat3.mat');
GV=load("GoodVars.mat");
GoodVars=GV.GoodVars;

%combine features from all sites
sitenums=[1:5]; %5 sites 
tempfeats=[feats{1};feats{2};feats{3};feats{4};feats{5}];
%set up featsB one field at a time
featsB.geom=[tempfeats(1).geom;tempfeats(2).geom;tempfeats(3).geom;tempfeats(4).geom;tempfeats(5).geom];
featsB.slope=[tempfeats(1).slope;tempfeats(2).slope;tempfeats(3).slope;tempfeats(4).slope;tempfeats(5).slope];
featsB.scan=[tempfeats(1).scan;tempfeats(2).scan;tempfeats(3).scan;tempfeats(4).scan;tempfeats(5).scan];
featsB.rough=[tempfeats(1).rough;tempfeats(2).rough;tempfeats(3).rough;tempfeats(4).rough;tempfeats(5).rough];
featsB.LABpoint=[tempfeats(1).LABpoint;tempfeats(2).LABpoint;tempfeats(3).LABpoint;tempfeats(4).LABpoint;tempfeats(5).LABpoint];
featsB.meanLAB=[tempfeats(1).meanLAB;tempfeats(2).meanLAB;tempfeats(3).meanLAB;tempfeats(4).meanLAB;tempfeats(5).meanLAB];
featsB.LABstd=[tempfeats(1).LABstd;tempfeats(2).LABstd;tempfeats(3).LABstd;tempfeats(4).LABstd;tempfeats(5).LABstd];
featsB.WaltonLAB=[tempfeats(1).WaltonLAB;tempfeats(2).WaltonLAB;tempfeats(3).WaltonLAB;tempfeats(4).WaltonLAB;tempfeats(5).WaltonLAB];
featsB.GLCM=[tempfeats(1).GLCM;tempfeats(2).GLCM;tempfeats(3).GLCM;tempfeats(4).GLCM;tempfeats(5).GLCM];
% featsB.LABtexture=[tempfeats(1).LABtexture;tempfeats(2).LABtexture;tempfeats(3).LABtexture;tempfeats(4).LABtexture]; %no Berretta features desired for testing 
featsB.intensity=[tempfeats(1).intensity;tempfeats(2).intensity;tempfeats(3).intensity;tempfeats(4).intensity;tempfeats(5).intensity];
featsB.truth=[tempfeats(1).truth;tempfeats(2).truth;tempfeats(3).truth;tempfeats(4).truth;tempfeats(5).truth];
featsB.points=[tempfeats(1).points;tempfeats(2).points;tempfeats(3).points;tempfeats(4).points;tempfeats(5).points];

%% Evaluate impact of training/testing ratio 
trainprop=[0.2:0.1:0.9];
for j=1:length(trainprop)
    for q=1:30
    numsamples=10000; 
    sinds=randsample(500000,numsamples);
    %set testing fraction: no overlap but classes not evenly sampled. 
    train1=sinds(1:floor(numsamples*trainprop(j)));
    test1=sinds(floor(numsamples*trainprop(j))+1:end);
    all = [featsB.geom,featsB.slope,featsB.intensity,featsB.meanLAB,featsB.LABstd,featsB.WaltonLAB,featsB.scan,featsB.rough,featsB.LABpoint,featsB.GLCM];
    %limit to variables from recursive feature elimination. 
    allt = featsB.truth;
    x_test=all(test1,GoodVars);
    x_train=all(train1,GoodVars);
    y_test = allt(test1,:);
    y_train = allt(train1,:);
    test_points=featsB.points(test1,:);
    rf = TreeBagger(30,x_train,y_train,'maxnumsplits',80);
    [y_pred,score] = predict(rf,x_test); %prediction
    y_prediction = str2double(y_pred); %accuracy evaluation
    stats1 = confusionmatStats(y_test,y_prediction);
    Fsc(q)=stats1.Fscore(1,1); %store F-Score
    %adjust decision threshold 
    new_pred=zeros(length(y_pred),1); %decision split threshold 
    threshold=0.65; %increase threshold to reduce false positives (false clasts) 
    new_pred(find(score(:,1)<threshold))=1; %large value in column 1 indicates more likely to be 0 class (clast)
    %get statistics - precision of class 0 should be over 0.9 for best results
    stats = confusionmatStats(y_test,new_pred);
    Prec(q)=stats.precision(1,1); % Store Precision
    end
Fsave{j}=Fsc;
Precsave{j}=Prec;
clear Prec Fsc
end

%% Reformat results for plotting and save data
trainprop=[0.2:0.1:0.9];
for j=1:length(trainprop)
    F2(:,j)=Fsave{j};
    P2(:,j)=Precsave{j};
end

save A1_RatioTest.mat F2 P2 trainprop
%% Plot 
% load('A1_RatioTest.mat'); %uncomment if desired
t=tiledlayout(1,2)
% ff=figure;
% ff.Position=[0 0 700 500];
nexttile
boxplot(P2(:,1:8),'Labels',trainprop)
ylabel('Precision at Decision Boundary=0.65','fontsize',14)
% xlabel('Training Proportion','fontsize',14)

nexttile
boxplot(F2(:,1:8),'Labels',trainprop)
ylabel('Fscore at Decision Boundary=0.5','fontsize',14)
xlabel(t,'Training Proportion','fontsize',14)

%% Save figure 
print('A1_RatioTest','-dpng', '-r400')