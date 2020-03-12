clear all; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

% plot raw signal for figure(1)
Un(:,:,:)=reshape(Undata(1,:),n,n,n);
close all, isosurface(X,Y,Z,abs(Un),0.4);
axis([-20 20 -20 20 -20 20]), grid on, drawnow;
xlabel('x-coordinates');
ylabel('y-coordinates');
zlabel('z-coordinates');
title('Raw Ultasound Signal');


%% 
% plot average signal at 5, 10, 20 for figure(2)
Uave = zeros(n,n,n);
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n);
Uave = Uave + fftn(Un);
end
Uave = Uave/20;
close all, isosurface(Kx,Ky,Kz,abs(Uave)/max(abs(Uave(:))),0.6);
% axis([-20 20 -20 20 -20 20]), grid on, drawnow;
xlabel('x-frequency');
ylabel('y-frequency');
zlabel('z-frequency');
title('Average Signal Frequency Domain (20 Signals)');


%% filter marble frequency
[max_val, position] = max(abs(Uave(:)));
[x, y, z] = ind2sub(size(Uave),position);
fx = Kx(x, y, z);
fy = Ky(x, y, z);
fz = Kz(x, y, z);
filter = exp(-2*pi*((Kx - fx).^2 + (Ky - fy).^2 + (Kz - fz).^2));

%%
target=zeros(3,20);
for j=1:20
    Un = reshape(Undata(j,:),n,n,n);
    Untf = fftn(Un).* filter;
    Unf = ifftn(Untf);
    [max_value, position] = max(abs(Unf(:)));
    [marblex, marbley, marblez] = ind2sub([n n n], position);
    target(1,j) = X(marblex, marbley, marblez);
    target(2,j) = Y(marblex, marbley, marblez);
    target(3,j) = Z(marblex, marbley, marblez);
end

figure()
plot3(target(1,:), tar2get(2,:), target(3,:),'b'), grid on;
hold on; 
plot3(target(1,20), target(2,20), target(3,20),'r-o','LineWidth',[25]);
text(target(1,20), target(2,20), target(3,20), '   \leftarrow (-6.09, 4.22, -6.09)');
xlabel('x-coordinate');
ylabel('y-coordinate');
zlabel('z-coordinate');
title('Target Path and Final Location (Red)');