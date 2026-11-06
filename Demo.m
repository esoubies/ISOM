%--------------------------------------------------------------------------
% Demo script of the two ISOM reconstruction methods proposed in [1].
% 
% [1] Fourier Analysis of Interference Scanning Optical Probe Microscopy.
%     IEEE Transactions on Computational Imaging, vol. 11, pp. 1206-1216, 2025.
%     Emmanuel Soubies and Wolfgang Bacsa. 
%
% 2025 - Emmanuel Soubies (emmanuel.soubies@cnrs.fr)
%--------------------------------------------------------------------------
clear;close all;rng(1);
run('./Install')

%% Physical parameters
% - Incident wave
Lb = 0.5;              % Wavelength (in um)
theta_deg = 60;        % Angle of the incident wave with z axis in degree
phi_deg = 0;           % Angle of the incident wave with x axis in degree

% - Model mismatch in the reconstruction model
dtheta=0/180*pi;
dphi=0/180*pi;
dzim=0; 

% - Camera
cam_side = 90*Lb;      % Length of camera side (in um)
zim = 0.5*Lb;          % Distance from image plane to sensor (in um) 
scale =0;              % Scale for lateral distance between object and camera plane (in the wavevector direction)
N = 512;

% - Noise
sig_n = 0.1;           % std of the gaussian noise will be nn*max(ynonoise(:)) [If array -> loop]

% - Optim
MaxIter = 500;               % Maximum number of iterations
cv_tol = 1e-8;               % Tolerance for stopping criteria
RegParam =5e-8*sig_n.^2;     % Regularization (tuning) parameter [If sig_n is an array, need to be of the same length]
RegType = 1;                 % Regularization type
nonneg=1;                    % boolean true if nonnegativity is activated for full model reconstruction
alt_iter = 1;                % Number of alternate iterations
OTFatt=0;                    % boolean true is OTF attenuation is activated

%% Generate data
% - Object
im = double(loadtiff('phantom.tif')); im=im/max(im(:));

% - Incident wave
theta = (theta_deg*pi)/180;        % Conversion in radian
phi = (phi_deg*pi)/180;            % Conversion in radian
k0 = 2*pi/Lb;                      % Wavenumber
ki = [sin(theta)*cos(phi),sin(theta)*sin(phi),-cos(theta)]*k0;    % Incident wavevector

% - Camera
lat_disp = (zim/tan(theta))*scale;       % Lateral distance between object and camera plane (in the wavevector direction)
cam_shift=ki/norm(ki);
cam_shift=lat_disp*cam_shift/norm(cam_shift);  
px_sz=cam_side/N;

% - Fwd Model
H = BuildFwd(cam_side, N, zim, ki, cam_shift, k0, 0);

% - Data
ynonoise=H*im;
y=ynonoise+randn(size(ynonoise))*sig_n;

% - Display object and data
figure
subplot(2,2,1); imagesc(ynonoise); axis image; title('Noiseless data');
subplot(2,2,2); imagesc(log(1+abs(fftshift(fft2(ynonoise))))); axis image; title('FFT noiseless data');
subplot(2,2,3); imagesc(y); axis image; title('Noisy data');
subplot(2,2,4); imagesc(log(1+abs(fftshift(fft2(y))))); axis image; title('FFT noisy data');
colormap(flipud(gray))

%% Reconstruction
% -- For FRC computation 
p.beta=1;
freq=linspace(0,1/(px_sz*2),N/2+1);

% -- Reconstruction Full model
disp('== Full model reconstruction ==')
xFull = ReconsFullModel(cam_side,Lb,theta+dtheta,phi+dphi,zim+dzim,lat_disp,RegType,RegParam,nonneg,MaxIter,alt_iter,cv_tol,y,0);
PSNR_full = OptSNR(im,xFull);
frc = FSC(im,xFull,p);
idx = find(frc.frc<frc.T_hbit,1,'first');
res_full =  1/freq(idx);

% -- Joint deconvolution
disp('== Joint deconvolution ==')
[xJointDeconv,xc,xs] = JointDeconv(cam_side,Lb,theta+dtheta,phi+dphi,zim+dzim,lat_disp,RegType,RegParam,MaxIter,alt_iter,cv_tol,OTFatt,y,0);
PSNR_joint = OptSNR(im,xJointDeconv);
frc = FSC(im,xJointDeconv,p);
idx = find(frc.frc<frc.T_hbit,1,'first');
res_joint = 1/freq(idx);

%% Displays
% -- Metrics
disp(['PSNR Full: ',num2str(PSNR_full)])
disp(['PSNR JointDeconv: ',num2str(PSNR_joint)])
disp(['Res Full: ',num2str(res_full)])
disp(['Res JointDeconv: ',num2str(res_joint)]);

% -- Reconstructions
figure
subplot(2,2,1); imagesc(xFull); axis image; title('Full Model');
subplot(2,2,3); imagesc(log(1+abs(fftshift(fft2(xFull))))); axis image; title('FFT Full Model');
subplot(2,2,2); imagesc(xJointDeconv); axis image; title('Joint Deconv');
subplot(2,2,4); imagesc(log(1+abs(fftshift(fft2(xJointDeconv))))); axis image; title('FFT Joint Deconv');
colormap(flipud(gray))


