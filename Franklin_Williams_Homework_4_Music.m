close all; clear all;

%% Test 1
%load TalkingHeadsSamples.mat; Artist1 = THsamples;
%load WidespreadPanicSamples.mat; Artist2 = WSPsamples;
%load MozartSamples.mat; Artist3 = Msamples;
%fignum = 1; 

%% Test 2
%load GratefulDeadSamples.mat; Artist1 = GDsamples;
%load DeadAndCompanySamples.mat; Artist2 = DCsamples;
%load AlmostDeadSamples.mat; Artist3 = JRADsamples;
%fignum = 2; 

%% Test 3

% Bluegrass
load BillyStringsSamples.mat; A1 = BSsamples;
load GreenskySamples.mat; A2 = GSBGsamples;
load DelMcCouryBandSamples.mat; A3 = DMBsamples;
Artist1 = [A1(:,1:33) A2(:, 1:33) A3(:, 1:34)];

% HipHop
load HipHopSamples.mat; Artist2 = HHsamples;

% Jam
load WidespreadPanicSamples.mat; A1 = WSPsamples;
load GratefulDeadSamples.mat; A2 = GDsamples;
load AlmostDeadSamples.mat; A3 = JRADsamples;
Artist3 = [A1(:,1:33) A2(:, 1:33) A3(:, 1:34)];

fignum = 3

%%
Fs = 44100;
rng(512)
ctrain = [ones(80,1); 2*ones(80,1); 3*ones(80,1)];
ttrain = [ones(100,1); 2*ones(100,1); 3*ones(100,1)];
data = [Artist1(:,1:100) Artist2(:,1:100) Artist3(:,1:100)];

L = wmaxlev(length(data), 'haar');
cmat = [];
for i = 1:300
[c,l]=wavedec(data(:,i), L, 'haar');
cmat = [cmat c];
end

[u,s,v]=svd(cmat,0);
nbError = []; svmError = []; treeError = [];
for i = 1:100
Artist1 = v(:,randperm(100));
Artist2 = v(:,randperm(100)+100);
Artist3 = v(:,randperm(100)+200);

% Train and Test Data
train = [Artist1(1:80,:); Artist2(101:180,:); Artist3(201:280,:)];
test = [Artist1(81:100, :); Artist2(181:200,:); Artist3(281:300,:)];

t = templateNaiveBayes();
prenb = fitcecoc(train,ctrain,'CrossVal','on','Learners',t);
Nerr = kfoldLoss(prenb,'LossFun','ClassifErr');
nbError = [nbError Nerr];

t = templateSVM('Standardize',true);
Mdl = fitcecoc(train,ctrain,'Learners',t);
CVMdl = crossval(Mdl);
Serr = kfoldLoss(CVMdl);
svmError = [svmError Serr];

Mdl = fitctree(train,ctrain,'CrossVal','on');
Terr = kfoldLoss(Mdl);
treeError = [treeError Terr];
end

%% Histogram
figure(fignum)
hist1= histogram(nbError,8,'Normalization','pdf','EdgeColor', 'blue', 'FaceColor',  'blue');
title('Classification Method Error'); ylabel('Count'); xlabel('Model Error');
hold on
hist2 = histogram(svmError,10, 'EdgeColor', 'red', 'FaceColor',  'red', 'FaceAlpha', 0.2);
hist3 = histogram(treeError,10, 'EdgeColor', 'green', 'FaceColor',  'green', 'FaceAlpha', 0.2);
legend();
hold off
saveas(figure(fignum), strcat('Test ', num2str(fignum),' Histograms.png'), 'png')
%% Descriptive Statistics Table
Minimum = [min(nbError); min(svmError); min(treeError)];
Average = [mean(nbError); mean(svmError); mean(treeError)];
Median = [median(nbError); median(svmError); median(treeError)];
Max = [max(nbError); max(svmError); max(treeError);];
Range = [range(nbError); range(svmError); range(treeError);];
Std = [std(nbError); std(svmError); std(treeError)];
Model = {'Naive Bayes'; 'Support Vector Machine'; 'Decision Tree'}
T = table(Model, Minimum, Average, Median, Max, Range, Std)
writetable(T, strcat('Test ', num2str(fignum), 'Descriptive Statistics.xlsx'))