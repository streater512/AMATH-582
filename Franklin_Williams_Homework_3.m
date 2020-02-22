
%% load data
close all; clear all;
for N=[1:3]
    for C=[1:4]
        load(strcat('cam', num2str(N), '_', num2str(C), '.mat'))
    end
end

% simple harmonic
[xv1_1, xi1_1, yv1_1, yi1_1] = pointTracker(vidFrames1_1);
[xv2_1, xi2_1, yv2_1, yi2_1] = pointTracker(vidFrames2_1);
[xv3_1, xi3_1, yv3_1, yi3_1] = pointTracker(vidFrames3_1);

% noisy simple harmonic
[xv1_2, xi1_2, yv1_2, yi1_2] = pointTracker(vidFrames1_2);
[xv2_2, xi2_2, yv2_2, yi2_2] = pointTracker(vidFrames2_2);
[xv3_2, xi3_2, yv3_2, yi3_2] = pointTracker(vidFrames3_2);

% horizontal displacement
[xv1_3, xi1_3, yv1_3, yi1_3] = pointTracker(vidFrames1_3);
[xv2_3, xi2_3, yv2_3, yi2_3] = pointTracker(vidFrames2_3);
[xv3_3, xi3_3, yv3_3, yi3_3] = pointTracker(vidFrames3_3);

% horizontal displacement and rotation
[xv1_4, xi1_4, yv1_4, yi1_4] = pointTracker(vidFrames1_4);
[xv2_4, xi2_4, yv2_4, yi2_4] = pointTracker(vidFrames2_4);
[xv3_4, xi3_4, yv3_4, yi3_4] = pointTracker(vidFrames3_4);


%% align data

% simple harmonic
xv1_1 = xv1_1(1:226); xi1_1 = xi1_1(1:226); yv1_1 = yv1_1(1:226); yi1_1 = yi1_1(1:226);
xv2_1 = xv2_1(1:226); xi2_1 = xi2_1(1:226); yv2_1 = yv2_1(1:226); yi2_1 = yi2_1(1:226);
xv3_1 = xv3_1(1:226); xi3_1 = xi3_1(1:226); yv3_1 = yv3_1(1:226); yi3_1 = yi3_1(1:226);

% noisy simple harmonic
xv1_2 = xv1_2(1:226); xi1_2 = xi1_2(1:226); yv1_2 = yv1_2(1:226); yi1_2 = yi1_2(1:226);
xv2_2 = xv2_2(1:226); xi2_2 = xi2_2(1:226); yv2_2 = yv2_2(1:226); yi2_2 = yi2_2(1:226);
xv3_2 = xv3_2(1:226); xi3_2 = xi3_2(1:226); yv3_2 = yv3_2(1:226); yi3_2 = yi3_2(1:226);

% horizontal displacement
xv1_3 = xv1_3(1:230); xi1_3 = xi1_3(1:230); yv1_3 = yv1_3(1:230); yi1_3 = yi1_3(1:230);
xv2_3 = xv2_3(1:230); xi2_3 = xi2_3(1:230); yv2_3 = yv2_3(1:230); yi2_3 = yi2_3(1:230);
xv3_3 = xv3_3(1:230); xi3_3 = xi3_3(1:230); yv3_3 = yv3_3(1:230); yi3_3 = yi3_3(1:230);

% horizontal displacement and rotation
xv1_4 = xv1_4(1:392); xi1_4 = xi1_4(1:392); yv1_4 = yv1_4(1:392); yi1_4 = yi1_4(1:392);
xv2_4 = xv2_4(1:392); xi2_4 = xi2_4(1:392); yv2_4 = yv2_4(1:392); yi2_4 = yi2_4(1:392);
xv3_4 = xv3_4(1:392); xi3_4 = xi3_4(1:392); yv3_4 = yv3_4(1:392); yi3_4 = yi3_4(1:392);

% X matrix tests 1:4
X1=[xi1_1; yi1_1; xi2_1; yi2_1; xi3_1; yi3_1];
X2=[xi1_2; yi1_2; xi2_2; yi2_2; xi3_2; yi3_2];
X3=[xi1_3; yi1_3; xi2_3; yi2_3; xi3_3; yi3_3];
X4=[xi1_4; yi1_4; xi2_4; yi2_4; xi3_4; yi3_4];


%%
% simple harmonic plots
figure(1); sgtitle('Simple Harmonic Motion Data');
subplot(3,2,1), plot(xi1_1); title('X Location Camera 1');
subplot(3,2,2), plot(yi1_1); title('Y Location Camera 1');

subplot(3,2,3), plot(xi2_1); title('X Location Camera 2');
subplot(3,2,4), plot(yi2_1); title('Y Location Camera 2');

subplot(3,2,5), plot(xi3_1); title('X Location Camera 3');
subplot(3,2,6), plot(yi3_1); title('Y Location Camera 3');

% noisy simple harmonic plots
figure(2); sgtitle('Noisy Simple Harmonic Motion Data');
subplot(3,2,1), plot(xi1_2); title('X Location Camera 1');
subplot(3,2,2), plot(yi1_2); title('Y Location Camera 1');

subplot(3,2,3), plot(xi2_2); title('X Location Camera 2');
subplot(3,2,4), plot(yi2_2); title('Y Location Camera 2');

subplot(3,2,5), plot(xi3_2); title('X Location Camera 3');
subplot(3,2,6), plot(yi3_2); title('Y Location Camera 3');

% horizontal displacement plots
figure(3); sgtitle('Horizontal Displacement Motion Data');
subplot(3,2,1), plot(xi1_3); title('X Location Camera 1');
subplot(3,2,2), plot(yi1_3); title('Y Location Camera 1');

subplot(3,2,3), plot(xi2_3); title('X Location Camera 2');
subplot(3,2,4), plot(yi2_3); title('Y Location Camera 2');

subplot(3,2,5), plot(xi3_3); title('X Location Camera 3');
subplot(3,2,6), plot(yi3_3); title('Y Location Camera 3');

% horizontal displacement and rotation plots
figure(4); sgtitle('Horizontal Displacement and Rotation Motion Data');
subplot(3,2,1), plot(xi1_4); title('X Location Camera 1');
subplot(3,2,2), plot(yi1_4); title('Y Location Camera 1');

subplot(3,2,3), plot(xi2_4); title('X Location Camera 2');
subplot(3,2,4), plot(yi2_4); title('Y Location Camera 2');

subplot(3,2,5), plot(xi3_4); title('X Location Camera 3');
subplot(3,2,6), plot(yi3_4); title('Y Location Camera 3');

%% PCA

[uS,sS,vS,lambdaS,sigmaS,YS] = producePCProj(X1); %simple harmonic data
[uN,sN,vN,lambdaN,sigmaN,YN] = producePCProj(X2); % noisy simple harmonics PCA
[uH,sH,vH,lambdaH,sigmaH,YH] = producePCProj(X3); % horizontal displacement PCA
[uHR,sHR,vHR,lambdaHR,sigmaHR,YHR] = producePCProj(X4); % horizontal displacement and rotation PCA

lam=0.06; [R1rN,R2rN]=inexact_alm_rpca(X2.',lam); % noisy simple harmonics RPCA
lam=0.1; [R1rH,R2rH]=inexact_alm_rpca(X3.',lam); % horizontal RPCA
lam=0.1; [R1rHR,R2rHR]=inexact_alm_rpca(X4.',lam); % horizontal rotation RPCA

[urN,srN,vrN]=svd(R1rN); sigmarN = diag(srN); % noisy simple harmonics RPCA
[urH,srH,vrH]=svd(R1rH); sigmarH = diag(srH); % horizontal RPCA
[urHR,srHR,vrHR]=svd(R1rHR); sigmarHR = diag(srHR); % horizontal rotation RPCA

% plot PCA
plotPCA(sigmaS, 5, uS, "Simple Harmonic Data", 25000);
plotPCA(sigmaN, 6, uN, "Simple Harmonic Noisy Data", 25000);
plotPCA(sigmaH, 7, uH, "Horizontal Displacement Data", 50000);
plotPCA(sigmaHR, 8, uHR, "Horizontal Displacement and Rotation Data", 8000);

plotPCA(sigmarN, 9, urN, "Robust Simple Harmonic Noisy Data", 8000);
plotPCA(sigmarH, 10, urH, "Robust Horizontal Displacement Data", 13000);
plotPCA(sigmarHR, 11, urHR, "Robust Horizontal Displacement and Rotation Data", 13000);

%% calculate energy

energyS = sigmaS/sum(sigmaS); % simple harmonic POD
energyN = sigmaN/sum(sigmaN); % noisey simple harmonic POD
energyH = sigmaH/sum(sigmaH); % horizontal displacement POD 
energyHR = sigmaHR/sum(sigmaHR); % horizontal displacement and rotation POD

energyrN = sigmarN/sum(sigmarN); % robust simple harmonic noisy POD
energyrH = sigmarH/sum(sigmarH); % robust horizontal displacement POD
energyrHR = sigmarHR/sum(sigmarHR); % robust horizontal displacement and rotation POD


figure(12);
bar([energyS, energyN, energyrN, energyH, energyHR]); %plot energy
title('Percentage of Variance Explained by Mode');
legend('simple', 'noisy simple','robust noisy simple' ,'horizontal', 'horizontal rotation');

figure(13);
bar([energyrN, energyrH, energyrHR]);
title('Percentage of Variance Explained by Mode');
legend('robust noisy simple', 'robust horizontal', 'robust horizontal rotation');

%% Simple Harmonic Comparison
figure(14)
subplot(3,1,1), plot(sigmaS,'ko','Linewidth',[1.5]);
set(gca,'Fontsize',[13],'Xtick',[0:6]);
axis([0 6 0 20000]);
subplot(3,1,2), plot(sigmaN,'ko','Linewidth',[1.5]);
set(gca,'Fontsize',[13],'Xtick',[0:6]);
axis([0 6 0 12000]);
subplot(3,1,3), plot(sigmarN,'ko','Linewidth',[1.5]);
set(gca,'Fontsize',[13],'Xtick',[0:6]);
axis([0 6 0 12000]);


%% functions
function [x_vt, x_it, y_vt, y_it] = pointTracker(video)
x_vt = [];
x_it = [];
y_vt = [];
y_it = [];
    for f=1:size(video,4)
        frame=rgb2gray(video(:,:,:,f));
        [M,I]=max(frame(:));
        [x_i, y_i]=ind2sub(size(frame), I);
        [x_v, y_v]=ind2sub(size(frame), M);
        x_vt=[x_vt x_v];
        x_it=[x_it x_i];
        y_vt=[y_vt y_v];
        y_it=[y_it y_i];
    end

end

function [u,s,v,lambda,sigma,Y] = producePCProj(X)
[m,n]=size(X);   %  compute data size
mn=mean(X,2); %  compute mean for each row
X=X-repmat(mn,1,n);  % subtract mean
Cx=(1/(n-1))*X*X'; % Cov matrix
[u,s,v]=svd(Cx);  % perform the SVD
lambda=diag(s).^2;  % produce diagonal variances
sigma=diag(s);
Y=u'*Cx;  % produce the principal components projection
end

function plotPCA(sigma, fignum, u, grouptitle, ylim)
x=[1:1:size(u,1)];
figure(fignum);
sgtitle(grouptitle);

subplot(2,1,1), plot(sigma,'ko','Linewidth',[1.5]);
set(gca,'Fontsize',[13],'Xtick',[0:6]);
axis([0 6 0 ylim]);
title('PCA');
energy2=sum(sigma(1:2))/sum(sigma);
leg2 = strcat('Sum Mode 1:2 energy=', num2str(energy2));
legend(leg2, 'Location', 'NorthEast') ;

subplot(2,1,2), semilogy(sigma,'ko','Linewidth',[1.5]);
axis([0 6 0 10^(5)]);
set(gca,'Fontsize',[13],'Ytick',[10^(-15) 10^(-10) 10^(-5) 10^0 10^5],...
      'Xtick',[0:6]);
title('Log PCA');
end