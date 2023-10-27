% Description: The following code applies the classifier+clustering 
% workflow to subregions of the front half of the debris flow with no
% spatial overlap between training and testing regions (cross validation)
%
% Outputs: Tiled layout with 5 DBloops subfigures comparing estimated and true clast sizes
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%%
clear 
load('Feat3.mat'); 
GV=load('GoodVars.mat');
GoodVars=GV.GoodVars;
numsites=5;

%read in point clouds for each subregion and get GSDs
for i=1:numsites
   fcloud=readmatrix(sprintf('ab%d.txt',i));%read in subregion
   ccloud=fcloud(~isnan(fcloud(:,7)),[1 2 3 7]);%clast points with clast index, non-clast points have NaN index
   gsdt{i} = GetGSD(ccloud,99.5); %get GSD from indexed clast cloud 
   pcloud=fcloud(:,1:3); %cloud of X,Y,Z point only 
   cloud{i}=pcloud; %save X Y Zclouds for later merging and interpolation
end

%initialize figure for tiledlayout of GSDs
f=figure;
f.Position = [0 0 900 900];
t=tiledlayout(3,2);

%% classify + cluster different combinations of 5 subregions (test on 4 train on 1)
for n=1:5
    featsA=feats{n};
    sitenums=1:5; %5 is just numsites, locked for testing code 
    Blist=setdiff(sitenums,n); %all site indexes NOT in training data to be classified and clustered
    tempfeats=[feats{Blist(1)};feats{Blist(2)};feats{Blist(3)};feats{Blist(4)}]; %combined feature set from indexes above-order doesn't matter at the moment 
    featsB.geom=[tempfeats(1).geom;tempfeats(2).geom;tempfeats(3).geom;tempfeats(4).geom]; %set up featsB one field at a time
    featsB.slope=[tempfeats(1).slope;tempfeats(2).slope;tempfeats(3).slope;tempfeats(4).slope];
    featsB.scan=[tempfeats(1).scan;tempfeats(2).scan;tempfeats(3).scan;tempfeats(4).scan];
    featsB.rough=[tempfeats(1).rough;tempfeats(2).rough;tempfeats(3).rough;tempfeats(4).rough];
    featsB.LABpoint=[tempfeats(1).LABpoint;tempfeats(2).LABpoint;tempfeats(3).LABpoint;tempfeats(4).LABpoint];
    featsB.meanLAB=[tempfeats(1).meanLAB;tempfeats(2).meanLAB;tempfeats(3).meanLAB;tempfeats(4).meanLAB];
    featsB.LABstd=[tempfeats(1).LABstd;tempfeats(2).LABstd;tempfeats(3).LABstd;tempfeats(4).LABstd];
    featsB.WaltonLAB=[tempfeats(1).WaltonLAB;tempfeats(2).WaltonLAB;tempfeats(3).WaltonLAB;tempfeats(4).WaltonLAB];
    featsB.GLCM=[tempfeats(1).GLCM;tempfeats(2).GLCM;tempfeats(3).GLCM;tempfeats(4).GLCM];
    featsB.LABtexture=[tempfeats(1).LABtexture;tempfeats(2).LABtexture;tempfeats(3).LABtexture;tempfeats(4).LABtexture];
    featsB.intensity=[tempfeats(1).intensity;tempfeats(2).intensity;tempfeats(3).intensity;tempfeats(4).intensity];
    featsB.truth=[tempfeats(1).truth;tempfeats(2).truth;tempfeats(3).truth;tempfeats(4).truth];
    featsB.points=[tempfeats(1).points;tempfeats(2).points;tempfeats(3).points;tempfeats(4).points];
    
    %Features ordered in groups of 50k after merging e.g. 50k clast/50k matrix/50k clast/50k matrix 
    %Desired preinterpolation density of 3000 testing points per square m
    %corresponds to 49800 test core points per region (testing data) 
    
    %set up even samples of clast/matrix for testing/training, randomly sampled
    c400=[1:50000 100001:150000 200001:250000 300001:350000];
    m400=[50001:100000 150001:200000 250001:300000 350001:400000];
    tst0=randsample(c400,49800*2); %even sampling of classes 
    tst1=randsample(m400,49800*2);
    testinds=[tst0 tst1];
   
    c100=[1:50000];
    m100=[50001:100000];
    trn0=randsample(c100,49800/2);
    trn1=randsample(m100,49800/2);
    traininds=[trn0 trn1];
    
    x_train = [featsA.geom,featsA.slope,featsA.scan,featsA.rough,featsA.intensity,featsA.LABpoint,featsA.meanLAB,featsA.LABstd,featsA.WaltonLAB,featsA.GLCM,featsA.LABtexture];
    x_test = [featsB.geom,featsB.slope,featsB.scan,featsB.rough,featsB.intensity,featsB.LABpoint,featsB.meanLAB,featsB.LABstd,featsB.WaltonLAB,featsB.GLCM,featsB.LABtexture];
    y_train = featsA.truth(traininds);
    y_test = featsB.truth(testinds);
    x_test = x_test(testinds,GoodVars);
    x_train = x_train(traininds,GoodVars);

    %create random forest model 
    rf = TreeBagger(30,x_train,y_train,'maxnumsplits',60);
    [y_pred,score] = predict(rf,x_test); %prediction
    y_prediction = str2double(y_pred); %accuracy evaluation
    new_pred=zeros(length(y_pred),1); %decision split threshold 
    threshold=0.65; %increase threshold to reduce false positives (false clasts) 
    new_pred(find(score(:,1)<threshold))=1;

    truepred = zeros(length(y_test),1);
    for i=1:length(truepred);
      if y_test(i)==0 && new_pred(i)==0 %true positive 
            truepred(i)=1;%if rock=0 then this is true rock
      elseif y_test(i)==0 && new_pred(i)==1 %false negative 
            truepred(i)=3;%false soil (+)
      elseif y_test(i)==1 && new_pred(i)==0%false positive 
            truepred(i)=4;%false rock (-)
      elseif y_test(i)==1 && new_pred(i)==1 %true negative
            truepred(i)=2; %true soil
      end
    end

    newcloud=[featsB.points(testinds,:) truepred]; %cloud to be interpolated
    FullCloud=[cloud{Blist(1)};cloud{Blist(2)};cloud{Blist(3)};cloud{Blist(4)}]; %full density cloud of testing region points only 
    interped=InterpLabel(FullCloud,newcloud);
    indtr =  find(interped(:,4) == 1); %true rock
    indfr =  find(interped(:,4) == 4); %false rock
    indr = [indtr; indfr]; %indexes of all rock in full cloud 
    clastcloud = interped(indr,1:3); 
    dens=1500;
    [clouda, idx] = DownCloud(clastcloud,dens); %downsampled cloud

    np=18; esfa=2.6; esfb=0.6; 
    [ea, eb] = EstimateEpsilon(np, esfa, esfb, clouda); %must run on downsampled cloud 

    gsdtrue=[gsdt{Blist(1)} gsdt{Blist(2)} gsdt{Blist(3)} gsdt{Blist(4)}]; %gsd for subregions 

    [gsd, cloudb, scsm, sc, prmse(i), frmse(i)]=DBloops(gsdtrue,clouda,np,ea,eb,2);
    titlestr=strcat("Train: ",num2str(n) ," Test: ",num2str(Blist(1)),",",num2str(Blist(2)),",",num2str(Blist(3)),",",num2str(Blist(4)));
    title(titlestr,'fontsize',14,'fontname','helvetica')
end

save F11_CrossValidation.mat prmse frmse %to later evaluate rmse if desired
%% formatting and legend
for i=1:5
af1=nexttile(i);
af1.FontSize=12;
af1.TitleFontSizeMultiplier=1.3;
xlim(af1,[0 1]);
end

ax = nexttile;
p_ax=ax.Position;
hold on
stairs(NaN,NaN,'LineWidth',2,'color','b')
stairs(NaN,NaN,'color','k','LineWidth',2)
plot(NaN,NaN,'color','g','LineWidth',2);
leg = legend({'True', 'Estimated','|Residuals|'});
p_leg=leg.Position;
delete(ax)
ax=axes('Position',[p_ax(1:2) 0 0]);
hold on
stairs(NaN,NaN,'LineWidth',2,'color','k') %true
stairs(NaN,NaN,'color','b','LineWidth',2) %estimated
plot(NaN,NaN,'color',[0.4660 0.6740 0.1880],'LineWidth',2); %residuals
leg = legend({'Training Data','DBloops Result','|Residuals|'});
leg.Location = 'none';
leg.FontSize = 18;
ax.Visible = false;
leg.Position=[0.63422222465608,0.158903176182147,0.219999995132287,0.140333336480459]; 
%final axes labels
t.XLabel.String='Grain Size (m)';
t.XLabel.FontSize=18;
t.YLabel.String='Proportion Passing';
t.YLabel.FontSize=18;
%% print figure 
print('F11_CrossValidation', '-dpng', '-r400')