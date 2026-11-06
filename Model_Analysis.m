%--------------------------------------------------------------------------
% Script that illustrate the Fourier analysis of ISOM made in [1].
% Specifically:
%   - Cell 1: set physical parameters.
%   - Cell 2: precomputations.
%   - Cell 3: Generate convolution kernels and modulation functions of the
%   model reformulation (10) in [1].
%   - Cell 4: Display these convolution kernels and modulation function
%   together with their Fourier transform.
%   - Cell 5: Computation of theoretical Fourier transforms of kernels
%   given by [1, Proposition 1]. Display comparisons with those obtained
%   numerically at Cell 4.
%   - Cell 6: Reproduce Figure 2 of [1].
%   - Cell 7: Reproduce Figure 3 of [1] which illustrates the equivalent.
%   model computation via first applying filters and then shift [1, Proposition 2].
%   - Cell 8: Reproduce Figure 4 of [1] which illustrates the evolution of
%   the derived theoretical resolution limit.
%
% [1] Fourier Analysis of Interference Scanning Optical Probe Microscopy.
%     IEEE Transactions on Computational Imaging, vol. 11, pp. 1206-1216, 2025.
%     Emmanuel Soubies and Wolfgang Bacsa. 
%
% 2025 - Emmanuel Soubies (emmanuel.soubies@cnrs.fr)
%--------------------------------------------------------------------------
clear; close all; clc;
run('./Install')

%% Cell 1: Physical parameters
% - Incident wave
Lb = 0.5;              % Wavelength (in um)
theta_deg = 60;        % Angle of the incident wave with z axis in degree
phi_deg = 0;           % Angle of the incident wave with x axis in degree

% - Camera
cam_side = 90*Lb;      % Length of camera side (in um)
h = 0.5*Lb;            % Distance from image plane to sensor (in um)
N = 512;               % Number of pixels 

% If one want to be at Shannon sampling rate (according the [1] with nu_eff=0)
% cam_side = Lb / (2*(1+sin((theta_deg*pi)/180))) * N; 

%% Cell 2: Precomputations
% - Physical parameters
theta = (theta_deg*pi)/180;                                       % Conversion in radian
phi = (phi_deg*pi)/180;                                           % Conversion in radian
k0 = 2*pi/Lb;                                                     % Wavenumber
ki = [sin(theta)*cos(phi),sin(theta)*sin(phi),-cos(theta)]*k0;    % Incident wavevector
px_sz=cam_side/N;                                                 % Pixel size

% - Apodisation mask (for display purposes)
[Xapo,Yapo]=meshgrid(linspace(-1,1,N),linspace(-1,1,N));
Apo=1./(1+exp(100*(abs(Xapo)-0.95)))./(1+exp(100*(abs(Yapo)-0.95)));

%% Cell 3: Generate convolution kernels (PSFs) and modulation functions
% - Enlarged grid for convolution kernels (to account for zero-padding)
rj=linspace(-cam_side,cam_side-cam_side/N,2*N);
[Xpsf,Ypsf]=meshgrid(rj);Ypsf=Ypsf(end:-1:1,:);
norm_r=sqrt(Xpsf.^2+Ypsf.^2+h^2);

% Sin convolution kernel (on enlarged grid)
hhs = sin((ki(1)*Xpsf + ki(2)*Ypsf) - k0*norm_r) ./ norm_r;
n_hs = norm(hhs(:),1);
hhs=hhs/n_hs;

% Cos convolution kernel (on enlarged grid) 
hhc = cos((ki(1)*Xpsf + ki(2)*Ypsf) - k0*norm_r) ./ norm_r;
n_hc = norm(hhc(:),1);
hhc=hhc/n_hc;

% Cos and Sin modulations functions (on original grid)
rrj=linspace(-cam_side/2,cam_side/2-cam_side/(2*N),N);
[X,Y]=meshgrid(rrj);Y=Y(end:-1:1,:);
ds= sin((ki(1)*X+ki(2)*Y + ki(3)*h));
dc = cos((ki(1)*X+ki(2)*Y + ki(3)*h));

%% Cell 4: Display convolution kernels and modulation functions
% -- Sampling points in Fourier space from -nu_s/2 to nu_s/2
fj=linspace(-1/(px_sz*2),1/(px_sz*2)-1/(px_sz*2*N),2*N); 
idx=round(linspace(1,2*N,5)); % subset of indexes for labeling axis
ffj=linspace(-1/(px_sz*2),1/(px_sz*2)-1/(px_sz*2*N),N);
jdx=round(linspace(1,N,5));   % subset of indexes for labeling axis

% -- Compute Fourier transforms of kernels
fft_hhs=fftshift(fft2(ifftshift(hhs)));
fft_hhc=fftshift(fft2(ifftshift(hhc)));
fft_ds=fftshift(fft2(Apo.*ds));
fft_dc=fftshift(fft2(Apo.*dc));

% -- Pre-computations
cen1=N+1+ki(1:2).*[1,-1]*px_sz/pi*N;   % Center 1
cen2=N+1-ki(1:2).*[1,-1]*px_sz/pi*N;   % Center 2
c_max=max(max(log(1+abs(fft_dc(:)))),max(log(1+abs(fft_ds(:)))));

% -- Display convolution kernels (real part, imaginary part and amplitude)
figure; colormap(flipud(gray))
% Sin kernel
subplot(2,4,1); imagesc(hhs); axis image; title('Sin kernel');colorbar;setAxis(idx,rj,1)
subplot(2,4,2); imagesc(real(fft_hhs)); axis image; title('Sin kernel (Real fft)');colorbar;setAxis(idx,fj,2);setAxisLabels(2*N);
subplot(2,4,3); imagesc(imag(fft_hhs)); axis image; title('Sin kernel (Imag fft)');colorbar;setAxis(idx,fj,2);setAxisLabels(2*N);
subplot(2,4,4); imagesc(abs(fft_hhs)); axis image; title('Sin kernel (Abs fft)');colorbar;setAxis(idx,fj,2);setAxisLabels(2*N);
hold all; plot(cen1(1),cen1(2),'or','markersize',5); plot(cen2(1),cen2(2),'or','markersize',5);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'b');
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen2(1)+xp,cen2(2)+yp,'b');


% Cos kernel
subplot(2,4,5); imagesc(hhc); axis image; title('Cos kernel');colorbar;setAxis(idx,rj,1)
subplot(2,4,6); imagesc(real(fft_hhc)); axis image; title('Cos kernel (Real fft)');colorbar;setAxis(idx,fj,2);setAxisLabels(2*N);
subplot(2,4,7); imagesc(imag(fft_hhc)); axis image; title('Cos kernel (Imag fft)');colorbar;setAxis(idx,fj,2);setAxisLabels(2*N);
subplot(2,4,8); imagesc(abs(fft_hhc)); axis image; title('Cos kernel (Abs fft)');colorbar;setAxis(idx,fj,2);setAxisLabels(2*N);
hold all; plot(cen1(1),cen1(2),'or','markersize',5); plot(cen2(1),cen2(2),'or','markersize',5);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'b');
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen2(1)+xp,cen2(2)+yp,'b');


% -- Display modulation patterns (real part, imaginary part and amplitude)
figure; 
% Sin modulation
subplot(2,4,1); imagesc(ds); axis image; title('Sin mod');colorbar;setAxis(jdx,rrj,1)
subplot(2,4,2); imagesc(log(1+abs(real(fft_ds)))); axis image; title('Sin mod (Abs Real fft)');colorbar;clim([0,c_max]);setAxis(jdx,ffj,2);setAxisLabels(N);
subplot(2,4,3); imagesc(log(1+abs(imag(fft_ds)))); axis image; title('Sin mod (Abs Imag fft)');colorbar;clim([0,c_max]);setAxis(jdx,ffj,2);setAxisLabels(N);
subplot(2,4,4); imagesc(log(1+abs(fft_ds))); axis image; title('Sin mod (Abs fft)');colorbar;clim([0,c_max]);setAxis(jdx,ffj,2);setAxisLabels(N);
hold all; plot(N/2+1+ki(1)/(2*pi)*px_sz*N,N/2+1-ki(2)/(2*pi)*px_sz*N,'or','markersize',5);
plot(N/2+1-ki(1)/(2*pi)*px_sz*N,N/2+1+ki(2)/(2*pi)*px_sz*N,'or','markersize',5);

% Cos modulation
subplot(2,4,5); imagesc(dc); axis image; title('Cos mod');colorbar;setAxis(jdx,rrj,1)
subplot(2,4,6); imagesc(log(1+abs(real(fft_dc)))); axis image; title('Cos mod (Abs Real fft)');colorbar;clim([0,c_max]);setAxis(jdx,ffj,2);setAxisLabels(N);
subplot(2,4,7); imagesc(log(1+abs(imag(fft_dc)))); axis image; title('Cos mod (Abs Imag fft)');colorbar;clim([0,c_max]);setAxis(jdx,ffj,2);setAxisLabels(N);
subplot(2,4,8); imagesc(log(1+abs(fft_dc))); axis image; title('Cos mod (Abs fft)');colorbar;clim([0,c_max]);setAxis(jdx,ffj,2);setAxisLabels(N);
hold all; plot(N/2+1+ki(1)/(2*pi)*px_sz*N,N/2+1-ki(2)/(2*pi)*px_sz*N,'or','markersize',5);
plot(N/2+1-ki(1)/(2*pi)*px_sz*N,N/2+1+ki(2)/(2*pi)*px_sz*N,'or','markersize',5);
colormap(flipud(gray))

%% Cell 5: Computation of theoretical Fourier transforms of kernels (see [1, Proposition 1]) and comparisons
% Preprocess
[Fj1,Fj2] = meshgrid(pi*fj*2,pi*fj*2); % fj defined in [m^-1] and we convert here to [rad/m] in order to apply our theoretical expression of the Fourier transform of the kernels.
norm_p = sqrt((ki(1) + Fj1).^2 + (ki(2) + Fj2).^2);
norm_m = sqrt((ki(1) - Fj1).^2 + (ki(2) - Fj2).^2);
gq = @(q,r)  (r < k0) .* 1i * (-1)^(q+1) * pi * sqrt(pi * h) ./ (sqrt(2) * (k0^2 - r.^2+eps).^(1/4)) .* besselh(-0.5,q,h * sqrt(k0^2 - r.^2+eps)) ...
           + (r > k0) .* pi ./ (sqrt(r.^2 - k0^2+eps)) .* exp(-h * sqrt(r.^2 - k0^2)) ;
% Above the +eps is to avoid NaN is for some pixels r = k0

% Theoretical Fourier transforms
fft_hhc_th = (gq(1,norm_p) + gq(2,norm_m))/(px_sz^2*n_hc); 
fft_hhs_th = 1i * (gq(1,norm_p) - gq(2,norm_m))/(px_sz^2*n_hs); 
id = abs(fft_hhc_th) > max(abs(fft_hhc(:))); fft_hhc_th(id) = fft_hhc_th(id) ./ abs(fft_hhc_th(id)) * max(abs(fft_hhc(:))); % Clip values too high (sampling points close to singularity, for the numerical FFT as the spacial kernel is cropped, FT is smooth)
id = abs(fft_hhs_th) > max(abs(fft_hhs(:))); fft_hhs_th(id) = fft_hhs_th(id) ./ abs(fft_hhs_th(id)) * max(abs(fft_hhs(:))); % Clip values too high (sampling points close to singularity, for the numerical FFT as the spacial kernel is cropped, FT is smooth)

% Comparisons Cos kernel
figure; 
s(1)=subplot(3,3,1); imagesc(real(fft_hhc_th));colorbar; title('Th Cos Kernel (Real fft)');clim([min(real(fft_hhc(:))),max(real(fft_hhc(:)))]);setAxis(idx,fj,2)
s(2)=subplot(3,3,2); imagesc(real(fft_hhc)); colorbar;title('Num Cos Kernel (Real fft)'); setAxis(idx,fj,2)
s(3)=subplot(3,3,3); imagesc(log(1+abs(real(fft_hhc_th) - real(fft_hhc)))); colorbar;title('Error (Real fft) - Log scale'); setAxis(idx,fj,2)
s(4)=subplot(3,3,4); imagesc(imag(fft_hhc_th));colorbar; title('Th Cos Kernel (Imag fft)');clim([min(imag(fft_hhc(:))),max(imag(fft_hhc(:)))]);setAxis(idx,fj,2)
s(5)=subplot(3,3,5); imagesc(imag(fft_hhc)); colorbar;title('Num Cos Kernel (Imag fft)'); setAxis(idx,fj,2)
s(6)=subplot(3,3,6); imagesc(log(1+abs(imag(fft_hhc_th)-imag(fft_hhc)))); colorbar;title('Error (Imag fft) - Log scale'); setAxis(idx,fj,2)
s(7)=subplot(3,3,7); imagesc(abs(fft_hhc_th));colorbar; title('Th Cos Kernel (Abs fft)');clim([min(abs(fft_hhc(:))),max(abs(fft_hhc(:)))]);setAxis(idx,fj,2)
s(8)=subplot(3,3,8); imagesc(abs(fft_hhc)); colorbar;title('Num Cos Kernel (Abs fft)'); setAxis(idx,fj,2)
colormap(flipud(gray));
s(9)=subplot(3,3,9); imagesc(log(1+abs(abs(fft_hhc_th)-abs(fft_hhc)))); colorbar;title('Error (Abs fft) - Log scale'); setAxis(idx,fj,2)
Link = linkprop(s,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);drawnow;

% Comparisons Sin kernel
figure; 
s(1)=subplot(3,3,1); imagesc(real(fft_hhs_th));colorbar; title('Th Sin Kernel (Real fft)');clim([min(real(fft_hhs(:))),max(real(fft_hhs(:)))]);setAxis(idx,fj,2)
s(2)=subplot(3,3,2); imagesc(real(fft_hhs)); colorbar;title('Num Sin Kernel (Real fft)'); setAxis(idx,fj,2)
s(3)=subplot(3,3,3); imagesc(log(1+abs(real(fft_hhs_th) - real(fft_hhs)))); colorbar;title('Error (Real fft) - Log scale'); setAxis(idx,fj,2)
s(4)=subplot(3,3,4); imagesc(imag(fft_hhs_th));colorbar; title('Th Sin Kernel (Imag fft)');clim([min(imag(fft_hhs(:))),max(imag(fft_hhs(:)))]);setAxis(idx,fj,2)
s(5)=subplot(3,3,5); imagesc(imag(fft_hhs)); colorbar;title('Num Sin Kernel (Imag fft)'); setAxis(idx,fj,2)
s(6)=subplot(3,3,6); imagesc(log(1+abs(imag(fft_hhs_th)-imag(fft_hhs)))); colorbar;title('Error (Imag fft) - Log scale'); setAxis(idx,fj,2)
s(7)=subplot(3,3,7); imagesc(abs(fft_hhs_th));colorbar; title('Th Sin Kernel (Abs fft)');clim([min(abs(fft_hhs(:))),max(abs(fft_hhs(:)))]);setAxis(idx,fj,2)
s(8)=subplot(3,3,8); imagesc(abs(fft_hhs)); colorbar;title('Num Sin Kernel (Abs fft)'); setAxis(idx,fj,2)
colormap(flipud(gray));
s(9)=subplot(3,3,9); imagesc(log(1+abs(abs(fft_hhs_th)-abs(fft_hhs)))); colorbar;title('Error (Abs fft) - Log scale'); setAxis(idx,fj/(2*pi),2)
Link = linkprop(s,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);drawnow;

%% Cell 6: Figure 2 of [1]
figure;
% Compute effective support
coef = 0.02;
neff = coef*abs(gq(1,0)+gq(2,2*norm(ki(1:2))));
veff = sqrt((Lambert_W(pi*h / (neff)) /(h))^2 + k0^2) - k0;
neff = neff / (px_sz^2*n_hc);

% Display Module Fourier of cos kernel
subplot(1,2,1); colormap(flipud(gray));
imagesc(log(1+log(1+log(1+log(1+log(1+abs(fft_hhc_th)))))));
title('$$|\hat{h}_c|$$','Interpreter','latex');clim([min(abs(fft_hhc(:))),max(abs(fft_hhc(:)))]);setAxis(idx,fj/(2*pi),2);
hold on; grid;
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); h1=plot(cen2(1)+xp,cen2(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'r--','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); h2=plot(cen2(1)+xp,cen2(2)+yp,'r--','Linewidth',2);
plot([0 2*N],[N N],'k--','Linewidth',1.8);
set(gca,'Fontsize',16);
hAxes = get(gca,'XAxis');
hAxes.TickLabelInterpreter = 'latex';
hAxes = get(gca,'YAxis');
hAxes.TickLabelInterpreter = 'latex';
set(gca,'YTick',[1 N 2*N]);set(gca,'YTickLabel',{'$$\frac{\nu_{s}}{2}$$', '0', '$$-\frac{\nu_{s}}{2}$$'});
set(gca,'XTick',[1 N 2*N]);set(gca,'XTickLabel',{'$$-\frac{\nu_{s}}{2}$$', '0', '$$\frac{\nu_{s}}{2}$$'});
h3=plot(cen1(1),cen1(2),'ok','markersize',6,'MarkerFaceColor','k'); plot(cen2(1),cen2(2),'ok','markersize',6,'MarkerFaceColor','k');
legend([h3,h1,h2],{'Centers $$\pm \mathbf{k}_r$$','Radius $$k_0$$','Radius $$k_0 + \nu_\mathrm{eff}$$'},'Interpreter','latex')

% Display profile
subplot(1,2,2);  hold on;
plot(fj,abs(fft_hhc_th(N+1,:)),'k','Linewidth',1.8);
plot(fj,ones(1,2*N)*neff,'k:','Linewidth',1.8);
plot([1 1]*norm(ki)*(1+sin(theta))/(2*pi),[1e-5 1e-1],'b-.','Linewidth',1.8);
plot(-[1 1]*norm(ki)*(1+sin(theta))/(2*pi),[1e-5 1e-1],'b-.','Linewidth',1.8);
plot([1 1]*(norm(ki)*(1+sin(theta))+veff)/(2*pi),[1e-5 1e-1],'r--','Linewidth',1.8);
plot(-[1 1]*(norm(ki)*(1+sin(theta))+veff)/(2*pi),[1e-5 1e-1],'r--.','Linewidth',1.8);
grid; set(gca,'Fontsize',16); axis([fj(1) fj(end) 1e-5 1e-1]);set(gca,'yscale','log');
set(gca,'YTickLabel',[]);
hAxes = get(gca,'XAxis');
hAxes.TickLabelInterpreter = 'latex';
set(gca,'XTick',[fj(1) 0 fj(end)]);set(gca,'XTickLabel',{'$$-\frac{\nu_{s}}{2}$$', '0', '$$\frac{\nu_{s}}{2}$$'});
xlabel('$$\nu_1$$','Interpreter','latex');
text(fj(1)*1.2,neff, '$$\eta_\mathrm{eff}$$','Interpreter','latex','Fontsize',16);
title('Line profile ($$\nu_2=0$$)','Interpreter','latex')

%% Cell 7: Figure 3 of [1]
% -- Application of the forward model on an image
% Read test image and compute its fft
im = double(loadtiff('phantom.tif'));
fftim=fftshift(fft2(im));

% Definition of cropping operator
sz = [2*N,2*N]; indxmin = [1,1]; indxmax = [N,N];
P = LinOpSelectorPatch(sz,indxmin,indxmax); 

% Modulations of image and ffts
dcim=dc.*im;
dsim=ds.*im;
fftdcim=fftshift(fft2((P'*dcim)));
fftdsim=fftshift(fft2((P'*dsim)));

% Convolution with kernels and ffts
fftcpart=fft_hhc_th.*fftdcim;
fftspart=fft_hhs_th.*fftdsim;
cpart=P*real(ifft2(ifftshift(fftcpart)));
spart=P*real(ifft2(ifftshift(fftspart)));

% Combination of cosine and sine parts
y=cpart-spart;
ffty=fftcpart-fftspart;

% Equivalent computation via first applying filters and then shift (see [1, Proposition 2])
norm_c = sqrt(Fj1.^2 + Fj2.^2);
fft_g1 = gq(1,norm_c)/(px_sz^2*n_hc); 
fft_g2 = gq(2,norm_c)/(px_sz^2*n_hs); 
mm = max(max(abs(fft_hhc(:))),max(abs(fft_hhs(:))));
id = abs(fft_g1) > mm; fft_g1(id) = fft_g1(id) ./ abs(fft_g1(id)) * mm; % Clip values too high (sampling points close to singularity, for the numerical FFT as the spacial kernel is cropped, FT is smooth)
id = abs(fft_g2) > mm; fft_g2(id) = fft_g2(id) ./ abs(fft_g2(id)) * mm; % Clip values too high (sampling points close to singularity, for the numerical FFT as the spacial kernel is cropped, FT is smooth)
fftg1part=imtranslate(fftshift(fft2((P'*im))).*fft_g1,-ki(1:2).*[1,-1]*px_sz/pi*N);
fftg2part=imtranslate(fftshift(fft2((P'*im))).*fft_g2,ki(1:2).*[1,-1]*px_sz/pi*N);
ffty_bis = exp(1i*ki(3)*h)*fftg2part + exp(-1i*ki(3)*h)*fftg1part;

% -- Generate Figure 3 of [1]
figure;
% Object
subplot(3,4,1); imagesc(im); axis image;hold on;  
title('$$\rho$$','Interpreter','latex'); set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'Fontsize',16);
subplot(3,4,5); 
imagesc(log(1+abs(fftim))); axis image; hold on;  
[xp,yp]=circle(norm(ki)*px_sz/pi*N/2*[1,1]); plot(N/2+1+xp,N/2+1+yp,'b-.','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N/2*[1,1]); plot(N/2+xp,N/2++yp,'r--','Linewidth',2);
arrow([0.4*N 1.05*N]/2-[ki(1)*px_sz/pi*N/2,0],[0.55*N 1.18*N]/2-[ki(1)*px_sz/pi*N/2,0],'Length',5,'Color',[0.4940 0.1840 0.5560]);
arrow([1.6*N 0.95*N]/2+[ki(1)*px_sz/pi*N/2,0],[1.45*N 0.82*N]/2+[ki(1)*px_sz/pi*N/2,0],'Length',5,'Color',[0.4940 0.1840 0.5560]);
title('$$|\hat \rho|$$','Interpreter','latex'); set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'Fontsize',16);

subplot(3,4,4); imagesc(y); axis image; setAxis(idx,rj,1); hold on; 
title('$$y$$','Interpreter','latex'); set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'Fontsize',16);

subplot(3,4,8); imagesc(log(1+abs(ffty))); axis image; setAxis(idx,fj/(2*pi),2); hold on; 
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); h1=plot(cen2(1)+xp,cen2(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'r--','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); h2=plot(cen2(1)+xp,cen2(2)+yp,'r--','Linewidth',2);
arrow([0.4*N 1.05*N],[0.55*N 1.18*N],'Length',5,'Color',[0.4940 0.1840 0.5560]);
arrow([1.6*N 0.95*N],[1.45*N 0.82*N],'Length',5,'Color',[0.4940 0.1840 0.5560]);
title('$$|\hat y|$$','Interpreter','latex'); set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'Fontsize',16);

% Cos Part
subplot(3,4,2); imagesc(cpart); axis image;  setAxis(idx,rj,1); hold on; 
title('$$h_c *(w_c \odot \rho)$$','Interpreter','latex');set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'Fontsize',16);

subplot(3,4,6); imagesc(log(1+abs(fftcpart))); axis image; setAxis(idx,fj,2);hold on; 
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); h1=plot(cen2(1)+xp,cen2(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'r--','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); h2=plot(cen2(1)+xp,cen2(2)+yp,'r--','Linewidth',2);
arrow([0.4*N 1.05*N],[0.55*N 1.18*N],'Length',5,'Color',[0.4940 0.1840 0.5560]);
arrow([1.6*N 0.95*N],[1.45*N 0.82*N],'Length',5,'Color',[0.4940 0.1840 0.5560]);
title('$$|\mathcal{F}\{h_c *(w_c \odot \rho)\}|$$','Interpreter','latex');set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'Fontsize',16);

% Sin Part
subplot(3,4,3); imagesc(spart); axis image;  setAxis(idx,rj,1); hold on; 
title('$$h_s *(w_s \odot \rho)$$','Interpreter','latex');set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'Fontsize',16);

subplot(3,4,7); imagesc(log(1+abs(fftspart))); axis image; setAxis(idx,fj,2);hold on; 
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); h1=plot(cen2(1)+xp,cen2(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'r--','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); h2=plot(cen2(1)+xp,cen2(2)+yp,'r--','Linewidth',2);
arrow([0.4*N 1.05*N],[0.55*N 1.18*N],'Length',5,'Color',[0.4940 0.1840 0.5560]);
arrow([1.6*N 0.95*N],[1.45*N 0.82*N],'Length',5,'Color',[0.4940 0.1840 0.5560]);
title('$$|\mathcal{F}\{h_s *(w_s \odot \rho)\}|$$','Interpreter','latex');set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
set(gca,'Fontsize',16);


subplot(3,4,10);
imagesc(log(1+abs(fftg1part))); axis image; setAxis(idx,fj/(2*pi),2); hold on; 
title('$$(\hat \rho \odot g_1(\|\cdot\|))(\cdot + \mathbf{k}_r)$$','Interpreter','latex'); set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'b-.','Linewidth',0.5);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen2(1)+xp,cen2(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'r--','Linewidth',0.5);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); plot(cen2(1)+xp,cen2(2)+yp,'r--','Linewidth',2);
arrow([1.6*N 0.95*N],[1.45*N 0.82*N],'Length',5,'Color',[0.4940 0.1840 0.5560]);
set(gca,'Fontsize',16);

subplot(3,4,11);
imagesc(log(1+abs(fftg2part))); axis image; setAxis(idx,fj/(2*pi),2); hold on; 
title('$$(\hat \rho \odot g_2(\|\cdot\|))(\cdot - \mathbf{k}_r)$$','Interpreter','latex'); set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'b-.','Linewidth',2);
[xp,yp]=circle(norm(ki)*px_sz/pi*N*[1,1]); plot(cen2(1)+xp,cen2(2)+yp,'b-.','Linewidth',0.5);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); plot(cen1(1)+xp,cen1(2)+yp,'r--','Linewidth',2);
[xp,yp]=circle((norm(ki)+veff)*px_sz/pi*N*[1,1]); plot(cen2(1)+xp,cen2(2)+yp,'r--','Linewidth',0.5);
arrow([0.4*N 1.05*N],[0.55*N 1.18*N],'Length',5,'Color',[0.4940 0.1840 0.5560]);
set(gca,'Fontsize',16);

colormap(flipud(gray))

%% Cell 8: Figure 4 of [1]: Evolution of theoretical resolution limit
neff_ref = 0.02*abs(gq(1,0)+gq(2,2*norm(ki(1:2)))); % this is the level considered in the previous figure.
v_eff = @(neff,zim) sqrt((Lambert_W(pi*zim / (neff)) ./(zim)).^2 + k0^2) - k0;
res = @(neff,zim) 1 ./ (1 + v_eff(neff,zim) /k0);
coef_list=[0.5  3  10];

figure;hold all;
zim_val = logspace(-1,0.5,50);

iii=1;
for c = coef_list
    res_val = res(c*neff_ref,zim_val*Lb);
    leg{iii}=['$$\eta_\mathrm{eff} =',num2str(c),'\,\eta_\mathrm{eff}^0$$'];
    plot(zim_val*Lb,res_val,'Linewidth',1.5);
    iii=iii+1;
end
legend(leg,'Interpreter','Latex','Location','southeast');
axis([0.17*Lb zim_val(end)*Lb 0.25 1]);
set(gca,'YTick',[0 0.25 0.5 0.75 1]);
hAxes = get(gca,'YAxis');
hAxes.TickLabelInterpreter = 'latex';
set(gca,'YTickLabel',{'0','$$\frac{\lambda}{4}$$','$$\frac{\lambda}{2}$$','$$\frac{3\lambda}{4}$$','$$\lambda$$'});
set(gca,'XTick',[0.17 0.5 1 2 3]*Lb);
hAxes = get(gca,'XAxis');
hAxes.TickLabelInterpreter = 'latex';
set(gca,'XTickLabel',{'$$\frac{\lambda}{6}$$','$$\frac{\lambda}{2}$$','$$\lambda$$','$$2 \lambda$$','$$3 \lambda$$'});
grid;xlabel('$$z_\mathrm{im}$$','Interpreter','Latex');ylabel('$$r_\mathrm{lim}$$','Interpreter','Latex');
set(gca,'Fontsize',14);


%% Utils functions
function setAxis(idx,f,domain)
if domain == 1
    set(gca,'xtick',idx,'xticklabel',(f(idx)));xlabel('$$r_1$$','Interpreter','Latex')
    set(gca,'ytick',idx,'yticklabel',(f(idx)));ylabel('$$r_2$$','Interpreter','Latex')
else
    set(gca,'xtick',idx,'xticklabel',(f(idx)));xlabel('$$\nu_1$$','Interpreter','Latex')
    set(gca,'ytick',idx,'yticklabel',(f(idx)));ylabel('$$\nu_2$$','Interpreter','Latex')
end
end

function [xp,yp]=circle(d)
ang=linspace(0,2*pi,100); 
xp=d(1)*cos(ang);yp=d(2)*sin(ang);
end

function setAxisLabels(N)
hAxes = get(gca,'XAxis');hAxes.TickLabelInterpreter = 'latex';
set(gca,'XTick',[1 N/2 N]);set(gca,'XTickLabel',{'$$-\frac{\nu_{s}}{2}$$', '0', '$$\frac{\nu_{s}}{2}$$'});
hAxes = get(gca,'YAxis');hAxes.TickLabelInterpreter = 'latex';
set(gca,'YTick',[1 N/2 N]);set(gca,'YTickLabel',{'$$\frac{\nu_{s}}{2}$$', '0', '$$-\frac{\nu_{s}}{2}$$'});
end