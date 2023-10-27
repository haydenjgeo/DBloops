% Description: The following code tests classification performance of a
% random forest model trained with selected features from RFE,
% training/testing on 10,000 total points 30x 
%
% Outputs: saves 1x30 vector of Fscores (clast class)named T1_FullWorkflow.mat
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%% random subsample of precalculated features
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


%%
%Sampling balanced by class but not position
q=[1 50001 100001 150001 200001 250001 300001 350001 400001 450001];
ind1=[]; ind0=[];
for i=1:length(q)
    if rem(i,2)~=0
    ind1s=[q(i):q(i)+49999];
    ind1=[ind1 ind1s]; %matrix
    else
    ind0s=[q(i):q(i)+49999];
    ind0=[ind0 ind0s]; %clast
    end
end

F=[]; %Fscore saver
for q=1:30

%testing fraction =0.8 (8000 points) from 250,000, so 4k each class
data0=randi(250000,1,5000); %even sampling of both classes (5k each)
y0=ind0(data0);
data1=randi(250000,1,5000);
y1=ind1(data1);
testi=[y0(1:4000) y1(1:4000)]; %no overlap in train/test
traini=[y0(4001:5000) y1(4001:5000)];

all = [featsB.geom,featsB.slope,featsB.intensity,featsB.meanLAB,featsB.LABstd,featsB.WaltonLAB,featsB.scan,featsB.rough,featsB.LABpoint,featsB.GLCM];

%limit to variables from recursive feature elimination 
allt = featsB.truth;
x_test=all(testi,GoodVars);
x_train=all(traini,GoodVars);
y_test = allt(testi,:);
y_train = allt(traini,:);

%create classifier and use on testing data
rf = TreeBagger(30,x_train,y_train,'maxnumsplits',60);
[y_pred,score] = predict(rf,x_test); %prediction
y_prediction = str2double(y_pred); %accuracy evaluation
stats1 = confusionmatStats(y_test,y_prediction);

F(q)=stats1.Fscore(1,1);
end
save T1_FullWorkflow.mat F