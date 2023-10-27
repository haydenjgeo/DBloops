% Description: The following code applies the classifier+clustering 
% workflow over the front half of the debris flow 
%
% Outputs: DBloops (1) and CompareClasts (2) figures comparing estimated and true clast sizes, merged into 3 panel layout 
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%% classification 
%Load precalculated features
load('Feat3.mat');
GV=load('GoodVars.mat');
GoodVars=GV.GoodVars;

%combine features from all sites
sitenums=[1:5]; %5 sites 
tempfeats=[feats{1};feats{2};feats{3};feats{4};feats{5}];
%set up featsB (testing data) one field at a time
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

%Downsampling - note core points are ordered as 50k clast/50k matrix/50k clast/50k matrix etc
numsamples=31125; %must be 50,000 or less - Feat3 contains 100k core points for 5 subregions, 50k matrix 50k clast
%note 31125 corresponds to 3000 ppm^2 testing data (3750 ppm^2 before 20/80 train/test split)

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
test_points=featsB.points(test1,:);


% Create classifier with precalculated features and classify "testing" data prior to clustering 
% don't worry about saving the indexes for your features - that's handled when you reinterpolate to the full density cloud 
rf = TreeBagger(30,x_train,y_train,'maxnumsplits',80);
[y_pred,score] = predict(rf,x_test); %prediction
y_prediction = str2double(y_pred); %accuracy evaluation
stats1 = confusionmatStats(y_test,y_prediction);

%result 7/26/23 - 81.58, 83.35 (Fscore, precision)

%adjust decision threshold for improved precision
new_pred=zeros(length(y_pred),1); %decision split threshold 
threshold=0.65; %increase threshold to reduce false positives (false clasts) 
new_pred(find(score(:,1)<threshold))=1; %large value in column 1 indicates likely 0 class (clast)
%get statistics - precision of class 0 should be over 0.9 for best results
stats = confusionmatStats(y_test,new_pred);

%result 7/26/23 - 74.07 Fscore, 92.16 precision after threshold adjustment. 

% assign labels to core points for extrapolating classifications to full density cloud 
% clasts are 0 (positive) and matrix material is 1 (negative) 
truepred = zeros(length(y_test),1);
for i=1:length(truepred)
  if y_test(i)==0 && new_pred(i)==0 %true positive 
        truepred(i)=1;%if rock=0 then this is true rock
  elseif y_test(i)==0 && new_pred(i)==1 %false negative 
        truepred(i)=3;%false soil 
  elseif y_test(i)==1 && new_pred(i)==0%false positive 
        truepred(i)=4;%false rock 
   elseif y_test(i)==1 && new_pred(i)==1 %true negative
        truepred(i)=2; %true soil
  end
end
newcloud=[test_points truepred]; %resulting cloud to be interpolated

%% clustering
fcloud=readmatrix('FH_NoUndersized_AllFields.txt'); %load in point cloud 
ccloud=fcloud(~isnan(fcloud(:,7)),[1 2 3 7]);%clast points with clast index, non-clast points have NaN index
gsdtrue = GetGSD(ccloud,99.5);
interped=InterpLabel(fcloud,newcloud);
%writematrix(interped,'interped.txt','Delimiter',','); % Uncomment to visualize extrapolated point cloud 
%Pretending we don't know what points are false positives to actually evaluate performance, so we include them 
indtr =  find(interped(:,4) == 1); %true rock
indfr =  find(interped(:,4) == 4); %false rock
indr = [indtr; indfr]; %indexes of all rock in full cloud 
%Cloud of just "clasts"
clastcloud = interped(indr,1:3); 
dens=1500;
%Downsample cloud for clustering, keep index of samples to compare to original cloud 
[clouda, idx] = DownCloud(clastcloud,dens);
%Define numpts, epsilon scaling factors
np=18; esfa=2.6; esfb=0.6; 
%Estimate epsilon base value
[ea, eb] = EstimateEpsilon(np, esfa, esfb, clouda); %must run on downsampled cloud 
%Run Clustering 
[gsd, cloudb, scsm, sc, prmse, frmse]=DBloops(gsdtrue,clouda,np,ea,eb,1);
gsdlength=length(gsd);

%results 7.26 full rmse 0.0317 partial rmse 0.0099

%Compare Clusters to Manually selected clasts 
%prepare "trucloud" at full resolution [x y z clastindex]
cloud=fcloud(:,[1 2 3 7]); %here my index 7 refers to the column of clastindexes 
%prepare cloud for compare clasts [x y z scsm sc pointindex]
cb6=indr(idx); %point index of samples
cloud1 = [cloudb cb6];
% clear i; %otherwise compareclasts3breaks 
rmse2=CompareClasts(cloud1, scsm, sc, cloud,1);

save F10_FullWorkflow.mat prmse frmse rmse2 gsdlength
%% combine figures produced by DBloops an CompareClasts into tiled layout
figlist=get(groot,'Children');
 
nf=figure;
nf.Position=[0 0 900 800];
tcl=tiledlayout(nf,'flow',"TileSpacing","compact")
figure(figlist(1))
ax=gca;
ax.Parent=tcl;
ax.Layout.Tile=2;
nexttile; 

figure(figlist(2))
ax=gca;
ax.Parent=tcl;
ax.Layout.Tile=3;

nexttile

figure(figlist(3))
ax=gca;
ax.Parent=tcl;
ax.Layout.Tile=1;
ax.Layout.TileSpan = [1 2];
nexttile

figure(4)
%%
annotation('textbox',[.16 .86 .04 .043],'String',' A','FitBoxToText','off','fontsize',14,'fontname','helvetica','fontweight','bold')
annotation('textbox',[.16 .41 .04 .043],'String',' B','FitBoxToText','off','fontsize',14,'fontname','helvetica','fontweight','bold')
annotation('textbox',[.58 .41 .04 .043],'String',' C','FitBoxToText','off','fontsize',14,'fontname','helvetica','fontweight','bold')
print('F10_FullWorkflow','-dpng', '-r400')
