close all; clear all;

load Cropped.mat
load UnCropped.mat

%% faces sample
figure(1);
 for j=1:9
 subplot(3,3,j);
 face1=uint8(reshape(CroppedMat(:,j), 192, 168));
 imshow(face1);
 end
 sgtitle('Example Cropped Face')
 saveas(figure(1),[pwd '/figures/figure(1).png'], 'png');

 
 figure(2);
 for j=1:9
 subplot(3,3,j);
 face2=uint8(reshape(UnCroppedMat(:,j), 243, 320));
 imshow(face2);
 end
 sgtitle('Example UnCropped Face')
 saveas(figure(2),[pwd '/figures/figure(2).png'], 'png');
 
  
%% Wavelet Representations
wavelets = {'haar', 'coif1', 'sym2', 'fk4', 'dmey'};
figure(3);
for w = 1:length(wavelets)
    [cod_cH1,cod_cV1,cod_edge] = wavelet_rep(CroppedMat,2,wavelets{w},192,168,w);
    subplot(5,4,w*4-3), imshow(uint8(cod_cH1));
    subplot(5,4,w*4-2), imshow(uint8(cod_cV1));
    subplot(5,4,w*4-1), imshow(uint8(cod_edge));
    subplot(5,4,w*4), imshow(uint8(reshape(CroppedMat(:,2), 192, 168)));
end
sgtitle('Wavelet Representation Examples Cropped Face');
saveas(figure(3),[pwd '/figures/figure(3).png'], 'png');

figure(4);
for w = 1:length(wavelets)
    [cod_cH1,cod_cV1,cod_edge] = wavelet_rep(UnCroppedMat,2,wavelets{w},243,320,w);
    subplot(5,4,w*4-3), imshow(uint8(cod_cH1));
    subplot(5,4,w*4-2), imshow(uint8(cod_cV1));
    subplot(5,4,w*4-1), imshow(uint8(cod_edge));
    subplot(5,4,w*4), imshow(uint8(reshape(UnCroppedMat(:,2), 243, 320)));
end
sgtitle('Wavelet Representation Examples UnCropped Face');
saveas(figure(4),[pwd '/figures/figure(4).png'], 'png');

%% Facial SVD
feature = 20;
figure(5)
for w = 1
    CroppedMat_wave = dc_wavelet(CroppedMat, wavelets{w}, 192, 168, 8064);
    [cU,cS,cV] = svd(CroppedMat_wave, 0);
    cU = cU(:,1:feature);
    
    subplot(2,2,w*4-3),ut1=reshape(cU(:,1),96,84);ut2=ut1(96:-1:1,:); pcolor(ut2);
    title('Principal Component 1');
    subplot(2,2,w*4-2),ut1=reshape(cU(:,2),96,84);ut2=ut1(96:-1:1,:); pcolor(ut2);
    title('Principal Component 2');
    subplot(2,2,w*4-1),ut1=reshape(cU(:,3),96,84);ut2=ut1(96:-1:1,:); pcolor(ut2);
    title('Principal Component 3');
    subplot(2,2,w*4),ut1=reshape(cU(:,4),96,84);ut2=ut1(96:-1:1,:); pcolor(ut2);
    title('Principal Component 4');
    set(gca,'Xtick',[],'Ytick',[]);
end
sgtitle('First Principal Components Cropped Faces')
saveas(figure(5),[pwd '/figures/figure(5).png'], 'png');

figure(6)
for w = 1
    UnCroppedMat_wave = dc_wavelet(UnCroppedMat, wavelets{w}, 243, 320, 19520);
    [uU,uS,uV] = svd(UnCroppedMat_wave, 0);
    uU = uU(:,1:feature);
    
    subplot(2,2,w*4-3),ut1=reshape(uU(:,1),122,160);ut2=ut1(122:-1:1,:); pcolor(ut2);
    title('Principal Component 1');
    subplot(2,2,w*4-2),ut1=reshape(uU(:,2),122,160);ut2=ut1(122:-1:1,:); pcolor(ut2);
    title('Principal Component 2');
    subplot(2,2,w*4-1),ut1=reshape(uU(:,3),122,160);ut2=ut1(122:-1:1,:); pcolor(ut2);
    title('Principal Component 3');
    subplot(2,2,w*4),ut1=reshape(uU(:,4),122,160);ut2=ut1(122:-1:1,:); pcolor(ut2);
    title('Principal Component 4');
    set(gca,'Xtick',[],'Ytick',[]);
end
sgtitle('First Principal Components UnCropped Faces');
saveas(figure(6),[pwd '/figures/figure(6).png'], 'png');

%%
cFeature = 1; uFeature = 1;
cExplained = 0; uExplained = 0;
cSigma = diag(cS); uSigma = diag(uS);
percent = [0.5, 0.75, 0.95];

cFeatureCount = []; uFeatureCount = [];

for p = 1:length(percent)
while cExplained < percent(p)
    cExplained = sum(cSigma(1:cFeature))/sum(cSigma);
    cFeature = cFeature+1;
end
cFeatureCount = [cFeatureCount, cFeature];

while uExplained < percent(p)
    uExplained = sum(uSigma(1:uFeature))/sum(uSigma);
    uFeature = uFeature+1;
end
uFeatureCount = [uFeatureCount, uFeature];
end
%% 
figure(7);
for i = 1:length(cFeatureCount)
    
subplot(3,2,i*2-1);
plot(cSigma,'ko','Linewidth',[1]);
sgtitle('Cropped Faces PCA');
title(strcat('Components Required to Explain  ', num2str(percent(i)*100), '%'), 'FontSize', [8])
hold on; plot(cSigma(1:cFeatureCount(i)), 'mo', 'Linewidth',[1]);
hold off; 
% set(gca,'Fontsize',[8],'Xlim',[0 length(cSigma)]);

subplot(3,2,i*2);
semilogy(cSigma,'ko','Linewidth',[1]); 
title(strcat('(Log Scale) Components Required to Explain  ', num2str(percent(i)*100), '%'), 'FontSize', [8])
hold on; semilogy(cSigma(1:cFeatureCount(i)),'mo','Linewidth',[1]);
hold off; 

% set(gca,'Fontsize',[8],'Xlim',[0 length(cSigma)]);
% sgtitle(strcat('Components Required to Explain  ', num2str(percent*100), '%')) 
end
saveas(figure(7),[pwd '/figures/figure(7).png'], 'png');

figure(8);
for i = 1:length(uFeatureCount)
    
subplot(3,2,i*2-1);
plot(uSigma,'ko','Linewidth',[1]);
sgtitle('UnCropped Faces PCA');
title(strcat('Components Required to Explain  ', num2str(percent(i)*100), '%'), 'FontSize', [8])
hold on; plot(uSigma(1:uFeatureCount(i)), 'mo', 'Linewidth',[1]);
hold off; 
% set(gca,'Fontsize',[8],'Xlim',[0 length(cSigma)]);

subplot(3,2,i*2);
semilogy(uSigma,'ko','Linewidth',[1]); 
title(strcat('(Log Scale) Components Required to Explain  ', num2str(percent(i)*100), '%'), 'FontSize', [8])
hold on; semilogy(uSigma(1:uFeatureCount(i)),'mo','Linewidth',[1]);
hold off; 

% set(gca,'Fontsize',[8],'Xlim',[0 length(cSigma)]);
% sgtitle(strcat('Components Required to Explain  ', num2str(percent*100), '%')) 
end
saveas(figure(8),[pwd '/figures/figure(8).png'], 'png');

%%
figure(9)
subplot(1, 2, 1); plot(100*[0:1/(length(cSigma)-1):1], 100*cumsum(cSigma)/sum(cSigma));
title('Cropped Faces'); ylabel('Explained (%)'); xlabel('Components Required (%)');
subplot(1, 2, 2); plot(100*[0:1/(length(uSigma)-1):1], 100*cumsum(uSigma)/sum(uSigma));
title('UnCropped Faces'); ylabel('Explained (%)'); xlabel('Components Required (%)');
saveas(figure(9),[pwd '/figures/figure(8).png'], 'png');

%%
figure(10)
for j=1:4
    subplot(4,2,2*j-1)
    plot(cV(1:150,j),'ko-')
    subplot(4,2,2*j)
    plot(uV(1:150,j),'ko-')
end
sgtitle('Projection onto First 3 POD Modes')
subplot(4,2,1); ylim([-0.1 0]); title('Cropped Faces');
subplot(4,2,2); ylim([-0.1 0]); title('UnCropped Faces');
saveas(figure(10),[pwd '/figures/figure(10).png'], 'png');

%% functions

function dcData = dc_wavelet(dcfile, w, x, y, nw)
[m,n]=size(dcfile); 
nbcol = size(colormap(gray),1);
 for i=1:n
     X=double(reshape(dcfile(:,i), x, y));
     [cA,cH,cV,cD]=dwt2(X,w);
     cod_cH1 = wcodemat(cH,nbcol);
     cod_cV1 = wcodemat(cV,nbcol);
     cod_edge=cod_cH1+cod_cV1;
     dcData(:,i)=reshape(cod_edge,nw,1);
 end
end

function [cod_cH1,cod_cV1,cod_edge] = wavelet_rep(mat,c,w,x,y,n)
X = double(reshape(mat(:,c), x, y));
[cA,cH,cV,cD]=dwt2(X,w);
nbcol = size(colormap(gray),1);
cod_cH1 = wcodemat(cH,nbcol);
cod_cV1 = wcodemat(cV,nbcol);
cod_edge=cod_cH1+cod_cV1;
end

