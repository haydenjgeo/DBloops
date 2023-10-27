% Description: This script compares Wolman pebble count results to the size distribution of all clasts manually trimmed
% from debris-flow point cloud data during creation of training data
%
% Outputs: figure comparing training data clast size distribution to Wolman count 
%
% Revision: R2022a
% Author: Hayden Jacobson, Colorado School of Mines/USGS
% Date: October 2023

%%
%get training data GSD
fcloud=readmatrix('FH_NoUndersized_AllFields.txt');
f=figure;
f.Position = [0 0 900 400];
ccloud=fcloud(~isnan(fcloud(:,7)),[1 2 3 7]);%clast points with clast index, non-clast points have NaN index
c2cloud=ccloud(:,[1 2 3]);%clast points no index
gsdtrue = GetGSD(ccloud,99.5);%GSD of manually derived training data
[f1,x1] = ecdf(gsdtrue);
stairs(x1,f1,'color',[33/256 49/256 77/256],'LineWidth',3);
hold on

%get Wolman Count GSD 
WC = readtable('WC42822.xlsx');
WC = WC(:,2);
WC = table2array(WC);
WC(isnan(WC)) = 2;
WC = WC/1000; %convert to meters passing (orig. mm)
WCu6ind = find(WC<=0.063);
WC(WCu6ind) = [];
[f2,x2]=ecdf(WC);
stairs(x2,f2,'color','r','LineWidth',3)

%formatting
settings.matlab.fonts = 'Tahoma'
set(gca,'yscale','linear','xscale','linear','xlim',[0 0.4]);
xlabel('Grain Diameter (m)','fontsize',16);
ylabel('Proportion Passing','fontsize',16);
legend('Training Data','Wolman Count','fontsize',14,'location','northeast');
xlim([0 max(gsdtrue)])

%%
print(gcf,'F8_Manual_WolmanCount.png','-dpng','-r500'); 
