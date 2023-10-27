% Description: The following code conducts recursive feature elimination by
% removing the most correlated, lowest performing feature until no
% correlation coefficients greater than 0.6 are observed between features
%
% Outputs: Feature indexes for retained features (vector of integer
% values) and figure displaying heatmap representing retained features. The
% F-Score of the clast class before and after RFE will appear in the
% command window.
%
% Note: The following features and scales are available in the precalculated feature set (Feat3.mat) 
%
% sample scales 1-9 [0.3 0.2 0.15 0.1 0.075 0.05 0.03 0.02 0.01]
%
% geom %10 features*SCALES; ie would expect 1 and 11 to be correlated 1-90
%     %geom = [omni,eigent,aniso,planar,linear,curv,scatt,vecs'];
% slope %4 features*SCALES 91-126
%    %[mean, std, skew, kurt] 
% intensity %4 features*SCALES 127-162
%    %[mean, std, skew, kurt]
% meanLAB %3 feature*SCALES 163-189
%    %[meanL meanA meanB]
% LABstd %3 feature*SCALES 190-216
%    %[stdL stdA stdB]
% WaltonLAB %9 feature*SCALES 219-297
%    %[mean-point point-min max-point].*[L A B]
% scan %1 feature*SCALES
%    %[inf entropy]
% rough %1 feature*SCALES
%    %[roughness]
% LABpoint %3 features, 1 scale (corepoints) 
%    %[L A B]
% GLCM %4 features, 1 scale
%    %Grey level co-occurence matrix [contrast correlation energy homogeneity] "Haralick features"
%LABtexture-192 features (Berretta features) are not assessed here, they don't improve the classifier
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%% load in precalculated features; all possible features at 9 scales 
clear
load('Feat3.mat')
%Combine all features (total 500k, 50k per class per subregion)
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

%% random subsample of precalculated feature points 
%even sampling with no overlap for training and testing, numsamples *5
numsamples=10000;
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
q2=[0:10]*numsamples+1; %ten is number of sites*number of classes
train1=[];
test1=[];
for u=1:length(q)
    test2=tind(q2(u):(q2(u)+numsamples*tsf-1));
    train2=tind((q2(u)+numsamples*tsf):(q2(u+1)-1));
    train1=[train1 train2];
    test1=[test1 test2];
end

all = [featsB.geom,featsB.slope,featsB.intensity,featsB.meanLAB,featsB.LABstd,featsB.WaltonLAB,featsB.scan,featsB.rough,featsB.LABpoint,featsB.GLCM];
% all = [featsB.geom,featsB.slope,featsB.intensity,featsB.meanLAB,featsB.LABstd,featsB.WaltonLAB,featsB.LABpoint,featsB.GLCM]; % contains no scan or rough 

allt = featsB.truth;
x_test=all(test1,:);
x_train=all(train1,:);
y_test = allt(test1,:);
y_train = allt(train1,:);

%% Recursive feature elimination (RFE)
A=x_train;
A(isnan(A))=0; %do this only once cause you want NaNs later, only changes a few values, but you get all correlations
val=1;
BadVarList=[];

%run iterative removal, one at a time based on correlation + importance 
while val>=0.6 %threshold for correlation coefficient above which variables will be cut 
%create model 
    rf = TreeBagger(30,A,y_train,'maxnumsplits',60,OOBpredictorimportance="on");
%find most correlated feature set used in model
    cor=corrcoef(A);%calc correlation coefficients for feature values
    cor=cor.*~(eye(size(cor)) & cor==1);%set diagonal to 0 w/ inv identify mat
    cor=abs(cor);%abs value for strong negative corrlelations 
    [val,cpos]=max(cor,[],'all'); %get max correlation and position 
    [cr,cc]=ind2sub(size(cor),cpos); %convert linear index to coordinate
%compare importance of two correlated variables, set less important to NaN
    if rf.OOBPermutedPredictorDeltaError(cr) > rf.OOBPermutedPredictorDeltaError(cc)
        A(:,cc)=NaN; BadVarList=[BadVarList cc];
    else A(:,cr)=NaN; BadVarList=[BadVarList cr];
    end
end

%% evalute performance of model w/ reduced feature set, print F-scores
%full feature set 
rf = TreeBagger(30,x_train,y_train,'maxnumsplits',80,OOBpredictorimportance="on");
B=x_test;
[y_pred,score] = predict(rf,B); %prediction
y_prediction = str2double(y_pred); %accuracy evaluation
stats = confusionmatStats(y_test,y_prediction);
stats.Fscore(1)

%RFE resultant feature set
GoodVars=(1:size(x_train,2));
GoodVars(BadVarList)=[];
A=x_train(:,GoodVars); %A is training features 
rf = TreeBagger(30,A,y_train,'maxnumsplits',30,OOBpredictorimportance="on");
B=x_test(:,GoodVars); %B is testing, with features matching A 
[y_pred,score] = predict(rf,B); %prediction
y_prediction = str2double(y_pred); %accuracy evaluation
stats1 = confusionmatStats(y_test,y_prediction);
stats1.Fscore(1)

save('Terpunkto\\Features\\GoodVars.mat','GoodVars');

%% plot heatmap of retained features after RFE only  
%34 multiscale features, 7 single scale features 
% mfs=32; only 32 when roughness and scan index removed (see line 82) 
load('GoodVars.mat')
mfs=34;
pfs=7;
%length of each feature group in order see line 44 (geom, slope, etc.)
grouplength=[10 4 4 3 3 9 1 1];
scales=9;

mlabels={'omnivariance','eigentropy','anisotropy','planarity','linearity','curvature','scatter','PCA1','PCA2','PCA3',...%multiscale feature labels
    'mean slope','std. dev. slope','skew slope','kurtosis slope',...
    'mean intensity','mean-point intensity','point-min. intensity','max.-point intensity',...
    'mean lightness','mean A','mean B','std. dev. lightness','std. dev. A','std. dev. B',...
    'mean-point lightness', 'point-min. lightness','max.-point lightness',...
    'mean-point A', 'point-min. A','max.-point A','mean-point B', 'point-min. B','max.-point B','scan index','roughness'};
plabels={'point lightness','point A','point B',... %point feature labels
    'GLCM contrast','GLCM correlation','GLCM energy','GLCM homogeneity'};
sc={'30', '20' ,'15', '10', '7.5' ,'5' ,'3' ,'2','1'}; %scales in cm

%build multiscale table
count=0; %counter that increases by one until total features (300+ reached)
for p=1:length(grouplength)
    fmat=[];
    for i=1:scales
        for j=1:grouplength(p)
            count=count+1;
            if ismember(count,GoodVars); 
                fmat(j,i)=1; 
            else; fmat(j,i)=0;
            end
        end
    end
    fmsave{p}=fmat;
end

%combine subtables to master table
fmat2=[];
for i=1:length(fmsave)
    fmat2=[fmat2; fmsave{i}];
end

%build corepoint feature table
fmatpoint=[];
for i=1:length(plabels)
    count=count+1;
    if ismember(count,GoodVars); fmatpoint(i,1)=1;
    else; fmatpoint(i,1)=0;
    end
end

ff=figure
ff.Position=[0 0 550 800]
tiledlayout(1,2)
nexttile

%multiscale feature heat map 
h=heatmap(fmat2);
Ylab=string(mlabels);
Xlab=string(sc);
h.YDisplayLabels=Ylab;
h.XDisplayLabels=Xlab;
h.Title='Multiscale Features Retained';
h.XLabel='Feature Scales (cm)';
h.ColorbarVisible='off';

%core point feature heat map
nexttile
h2=heatmap(fmatpoint);
Ylab2=string(plabels);
h2.YDisplayLabels=Ylab2;
h2.XDisplayLabels='Core Points (No Scales)';
h2.Title='Core Point Features Retained ';
h2.ColorbarVisible='off';

print('A4_RFE', '-dpng', '-r400')
