clear all; close all; clc

figurepath='/Users/Streater/Desktop/CFRM/AMATH 582/homework 2/figures/'

load handel
subplot(3,1,1);
v = y'/2;
v = v(1:length(v) -1);
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Hallelujah Chorus');

L=length(v)/Fs; n=length(v);
t2=linspace(0,L,n); t=t2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);

%% Gabor Filters Span (1)

figure(1)
vgt_spec=[];
a=[50, 25, 12.5, 6.25, 3.125, 0];
tslide=0:0.1:L
for i=1:length(a)
    for j=1:length(tslide)
        g=exp(-a(i)*(t-tslide(j)).^2); % Gabor
        vg=g.*v; vgt=fft(vg);
        vgt_spec=[vgt_spec; abs(fftshift(vgt))];    
    end
    subplot(3,2,i);
    sgtitle('Gabor Filters')
    pcolor(tslide, ks, vgt_spec.');, shading interp;
    set(gca,'Fontsize',[8]);
    colormap(hot);
    xlabel('time (t)');
    ylabel('frequency (\omega)');
    title(append('\alpha = ', num2str(a(i)), '; step = 0.1s'));
    vgt_spec=[];
end
saveas(figure(1), strcat(figurepath, 'figure1.png'), 'png');

%% Gabor Filters Step (2)

figure(2)
vgt_spec=[];
a=25;
step=[0.025, 0.05, 0.1, 0.25, 0.5, 0.75];
tslide={0:step(1):L;
    0:step(2):L;
    0:step(3):L;
    0:step(4):L;
    0:step(5):L;
    0:step(6):L};

for i=1:size(tslide, 1)
    for j=1:length(tslide{i})
        g=exp(-a*(t-tslide{i}(j)).^2); % Gabor
        vg=g.*v; vgt=fft(vg);
        vgt_spec=[vgt_spec; abs(fftshift(vgt))];    
    end
    subplot(3,2,i);
    sgtitle('Gabor Filters');
    pcolor(tslide{i}, ks, vgt_spec.'), shading interp;
    set(gca,'Fontsize',[8]);
    colormap(hot);
    xlabel('time (t)')
    ylabel('frequency (\omega)')
    title(append('\alpha = ', num2str(a), '; step = ', num2str(step(i)), 's'));
    vgt_spec=[];
end
saveas(figure(2), strcat(figurepath, 'figure2.png'), 'png')
%% Gaussian Span (3)
figure(3)
vgt_spec=[];
tslide=0:0.1:L
a=[50, 25, 12.5, 6.25, 3.125, 0]
for i=1:length(a)
    for j=1:length(tslide)
        % g=exp(-a(i)*(t-tslide(j)).^2); % Gabor
        % g=(1-(t-tslide(j)).^2).*exp(-a(i)*(t-tslide(j)).^2/2)/sqrt(2*pi) %Mexican Hat
        g=1/sqrt((a(i)+10^-6)/pi)*exp(-a(i)*(t-tslide(j)).^2); % Gaussian
        vg=g.*v; vgt=fft(vg);
        vgt_spec=[vgt_spec; abs(fftshift(vgt))];    
    end
    subplot(3,2,i);
    sgtitle('Gaussian Filters')
    pcolor(tslide, ks, vgt_spec.'), shading interp;
    set(gca,'Fontsize',[8]);
    colormap(hot);
    xlabel('time (t)')
    ylabel('frequency (\omega)')
    title(append('\alpha = ', num2str(a(i)), '; step = 0.1s'));
    vgt_spec=[];  
end
saveas(figure(3), strcat(figurepath, 'figure3.png'), 'png');

%% Gaussian Step (4)

figure(4)
vgt_spec=[];
a=12.5;
step=[0.025, 0.05, 0.1, 0.25, 0.5, 0.75];
tslide={0:step(1):L;
    0:step(2):L;
    0:step(3):L;
    0:step(4):L;
    0:step(5):L;
    0:step(6):L};

for i=1:size(tslide, 1)
    for j=1:length(tslide{i})
        % g=exp(-a*(t-tslide{i}(j)).^2); % Gabor
        % g=(1-(t-tslide{i}(j)).^2).*exp(-a*(t-tslide{i}(j)).^2/2)/sqrt(2*pi) %Mexican Hat
        g=1/sqrt(a/pi)*exp(-a*(t-tslide{i}(j)).^2); % Gaussian
        vg=g.*v; vgt=fft(vg);
        vgt_spec=[vgt_spec; abs(fftshift(vgt))];    
    end
    subplot(3,2,i);
    sgtitle('Gaussian Filters');
    pcolor(tslide{i}, ks, vgt_spec.'), shading interp;
    set(gca,'Fontsize',[8]);
    colormap(hot);
    xlabel('time (t)')
    ylabel('frequency (\omega)')
    title(append('\alpha = ', num2str(a), '; step = ', num2str(step(i)), 's'));
    vgt_spec=[];
end
saveas(figure(4), strcat(figurepath, 'figure4.png'), 'png');

%% Mexican Hat Span (5)

figure(5)
vgt_spec=[];
tslide=0:0.1:L
a=[50, 25, 12.5, 6.25, 3.125, 0]
for i=1:length(a)
    for j=1:length(tslide)
        % g=exp(-a(i)*(t-tslide(j)).^2); % Gabor
        g=(1-(t-tslide(j)).^2).*exp(-a(i)*(t-tslide(j)).^2/2)/sqrt(2*pi); %Mexican Hat
        vg=g.*v; vgt=fft(vg);
        vgt_spec=[vgt_spec; abs(fftshift(vgt))];    
    end
    subplot(3,2,i);
    sgtitle('Mexican Hat Filters')
    pcolor(tslide, ks, vgt_spec.'), shading interp;
    set(gca,'Fontsize',[8]);
    colormap(hot);
    xlabel('time (t)')
    ylabel('frequency (\omega)')
    title(append('\alpha = ', num2str(a(i)), '; step = 0.1s'));
    vgt_spec=[];  
end
saveas(figure(5), strcat(figurepath, 'figure5.png'), 'png');

%% Mexican Hat Filter Step (6)

figure(6)
vgt_spec=[];
a=12.5;
step=[0.025, 0.05, 0.1, 0.25, 0.5, 0.75];
tslide={0:step(1):L;
    0:step(2):L;
    0:step(3):L;
    0:step(4):L;
    0:step(5):L;
    0:step(6):L};

for i=1:size(tslide, 1)
    for j=1:length(tslide{i})
        % g=exp(-a*(t-tslide{i}(j)).^2); % Gabor
        g=(1-(t-tslide{i}(j)).^2).*exp(-a*(t-tslide{i}(j)).^2/2)/sqrt(2*pi); %Mexican Hat
        vg=g.*v; vgt=fft(vg);
        vgt_spec=[vgt_spec; abs(fftshift(vgt))];    
    end
    subplot(3,2,i);
    sgtitle('Mexican Hat Filters');
    pcolor(tslide{i}, ks, vgt_spec.'), shading interp;
    set(gca,'Fontsize',[8]);
    colormap(hot);
    xlabel('time (t)')
    ylabel('frequency (\omega)')
    title(append('\alpha = ', num2str(a), '; step = ', num2str(step(i)), 's'));
    vgt_spec=[];
end
saveas(figure(6), strcat(figurepath, 'figure6.png'), 'png');

%% Haar Filter
figure(7)
vgt_spec=[];
tslide=0:0.1:L;
a=[50, 25, 12.5, 6.25, 3.125, 0];
for i=1:length(a)
    for j=1:length(tslide)
        % g=exp(-a(i)*(t-tslide(j)).^2); % Gabor
        % g=(1-(t-tslide(j)).^2).*exp(-a(i)*(t-tslide(j)).^2/2)/sqrt(2*pi); %Mexican Hat
        g=heaviside(t-tslide(j))-2*heaviside(t-(tslide(j)+((a(i)+10^-6)/4)))+heaviside(t-(tslide(j)+2*((a(i)+10^-6)/4))); %Shannon
        vg=g.*v; vgt=fft(vg);
        vgt_spec=[vgt_spec; abs(fftshift(vgt))];    
    end
    subplot(3,2,i);
    sgtitle('Haar Filters')
    pcolor(tslide, ks, vgt_spec.'), shading interp;
    set(gca,'Fontsize',[8]);
    colormap(hot);
    xlabel('time (t)')
    ylabel('frequency (\omega)')
    title(append('\alpha = ', num2str(a(i)), '; step = 0.1s'));
    vgt_spec=[];  
end
saveas(figure(7), strcat(figurepath, 'figure7.png'), 'png');

%%
figure(8);
vgt_spec=[];
a=12.5;
step=[0.025, 0.05, 0.1, 0.25, 0.5, 0.75];
tslide={0:step(1):L;
    0:step(2):L;
    0:step(3):L;
    0:step(4):L;
    0:step(5):L;
    0:step(6):L};

for i=1:size(tslide, 1)
    for j=1:length(tslide{i})
        % g=exp(-a*(t-tslide{i}(j)).^2); % Gabor
        % g=(1-(t-tslide{i}(j)).^2).*exp(-a*(t-tslide{i}(j)).^2/2)/sqrt(2*pi); %Mexican Hat
        g=heaviside(t-tslide{i}(j))-2*heaviside(t-(tslide{i}(j)+((a+10^-6)/4)))+heaviside(t-(tslide{i}(j)+2*((a+10^-6)/4))); %Shannon
        vg=g.*v; vgt=fft(vg);
        vgt_spec=[vgt_spec; abs(fftshift(vgt))];    
    end
    subplot(3,2,i);
    sgtitle('Haar  Filters');
    pcolor(tslide{i}, ks, vgt_spec.'), shading interp;
    set(gca,'Fontsize',[8]);
    colormap(hot);
    xlabel('time (t)')
    ylabel('frequency (\omega)')
    title(append('\alpha = ', num2str(a), '; step = ', num2str(step(i)), 's'));
    vgt_spec=[];
end
saveas(figure(8), strcat(figurepath, 'figure8.png'), 'png');
%% Piano
figure(999999)
subplot(2,1,1)
tr_piano=16; % record time in seconds
v=audioread('music1.wav').';
Fs=length(v)/tr_piano;
plot((1:length(v))/Fs,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow

figure(8)
pong = bandpass(v, [200 500], Fs);
bandpass(v, [200 500], Fs)
figure(9)
pspectrum(pong,Fs,'spectrogram', 'FrequencyLimits',[1 500],'FrequencyResolution',5,'MinThreshold',-35)

%% Recorder 
tr_rec=14; % record time in seconds
v=audioread('music2.wav').'; Fs=length(v)/tr_rec;
plot((1:length(v))/Fs,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
% p8 = audioplayer(v,Fs); playblocking(p8);

figure(10)
pong = bandpass(v, [2000 3000], Fs);
bandpass(v, [2000 3000], Fs)
figure(11)
pspectrum(pong,Fs,'spectrogram', 'FrequencyLimits',[700 1200],'FrequencyResolution',15,'MinThreshold',-125)
