%% This section loops through 5 subsections with our "Generalization" approach 100k training 400k testing 
clear 
load('Feat3.mat'); 
GV=load('GoodVars.mat');
GoodVars=GV.GoodVars;
numsites=5;

%%
F=[];

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

    q=[1 50001 100001 150001 200001 250001 300001 350001];
    ind1=[]; ind0=[];
    for i=1:length(q)
        if rem(i,2)~=0
            ind1s=[q(i):q(i)+49999];
            ind1=[ind1 ind1s];
        else
            ind0s=[q(i):q(i)+49999];
            ind0=[ind0 ind0s];
        end
    end
    
    for m=1:30
        data0=randi(200000,1,4000); %testing data - clast
        y0=ind0(data0);
        data1=randi(200000,1,4000); %testing data - matrix
        y1=ind1(data1);
        %5k so no overlap in train/test
        testinds=[y0 y1];
        traininds=[randi([1 50000],1,1000) randi([50001 100000],1,1000)]; % even samples of training data

        x_train = [featsA.geom,featsA.slope,featsA.scan,featsA.rough,featsA.intensity,featsA.LABpoint,featsA.meanLAB,featsA.LABstd,featsA.WaltonLAB,featsA.GLCM,featsA.LABtexture];
        x_test = [featsB.geom,featsB.slope,featsB.scan,featsB.rough,featsB.intensity,featsB.LABpoint,featsB.meanLAB,featsB.LABstd,featsB.WaltonLAB,featsB.GLCM,featsB.LABtexture];
        y_train = featsA.truth(traininds);
        y_test = featsB.truth(testinds);
        x_test = x_test(testinds,GoodVars);
        x_train = x_train(traininds,GoodVars);


        rf = TreeBagger(30,x_train,y_train,'maxnumsplits',60);
        [y_pred,score] = predict(rf,x_test); %prediction
        y_prediction = str2double(y_pred); %accuracy evaluation
        stats1 = confusionmatStats(y_test,y_prediction);
        F(n,m)=stats1.Fscore(1,1);
    end
end

save T1_CrossValidation.mat F

