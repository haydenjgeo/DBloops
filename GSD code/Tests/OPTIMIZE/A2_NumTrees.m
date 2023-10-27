% Description: The following code compares the out-of-bag permuted
% predictor delta error in a Treebagger model at different numbers of
% leaves, only requiring the creation of one classifier
%
% Outputs: Plots showing relationship between OOBPPDE and number of leaves
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%% random subsample of precalculated feature points
%this one runs pretty quick. One treebagger run. 
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

%even sampling with no overlap for training and testing, numsamples *5
numsamples=10000; %must be 50,000 or less - Feat3 contains 100k core points for 5 subregions, 50k matrix 50k clast

%vector of start point of each class in each subregion
q=[1 50001 100001 150001 200001 250001 300001 350001 400001 450001];
tind=[];
for i=1:length(q)
    tinds=randi([q(i) q(i)+49999],1,numsamples);
    tind=[tind tinds];
end
%you now have a vector of subsampled features 10*numsamples long

%set testing fraction 
tsf=0.8;

%extract 20% of data for testing, leave the rest 
q2=[0:10]*numsamples+1;
train1=[];
test1=[];
for u=1:length(q)
    test2=tind(q2(u):(q2(u)+numsamples*tsf-1));
    train2=tind((q2(u)+numsamples*tsf):(q2(u+1)-1));
    train1=[train1 train2];
    test1=[test1 test2];
end

all = [featsB.geom,featsB.slope,featsB.intensity,featsB.meanLAB,featsB.LABstd,featsB.WaltonLAB,featsB.scan,featsB.rough,featsB.LABpoint,featsB.GLCM];

%limit to variables from recursive feature elimination 
allt = featsB.truth;
x_test=all(test1,GoodVars);
x_train=all(train1,GoodVars);
y_test = allt(test1,:);
y_train = allt(train1,:);


%% Create classifier from precalculated features and classify "testing" data prior to clustering 
rf = TreeBagger(100,x_train,y_train,'maxnumsplits',80,OOBpredictorimportance="on");

%% Create figure 
ff=figure;
ff.Position=[0 0 700 500];
plot(oobError(rf))
xlabel("Number of Grown Trees")
ylabel("Out-of-Bag Classification Error")

%% Save figure
print('A2_OOBvsNumTrees','-dpng', '-r400')

