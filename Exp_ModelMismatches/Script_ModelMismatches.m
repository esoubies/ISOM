%--------------------------------------------------------------------------
% Script that reproduce the robustness to model mismatches
% presented in Figure 8 of [1].
% 
% [1] Fourier Analysis of Interference Scanning Optical Probe Microscopy.
%     IEEE Transactions on Computational Imaging, vol. 11, pp. 1206-1216, 2025.
%     Emmanuel Soubies and Wolfgang Bacsa. 
%
% 2025 - Emmanuel Soubies (emmanuel.soubies@cnrs.fr)
%--------------------------------------------------------------------------
clear;close all;rng(1);
run('../Install')

%% Parameters
% - Incident wave
Lb = 0.5;              % Wavelength (in um)
theta_deg = 60;        % Angle of the incident wave with z axis in degree
phi_deg = 0;           % Angle of the incident wave with x axis in degree

% - Camera
cam_side = 90*Lb;      % Length of camera side (in um)
zim = 0.5*Lb;          % Distance from image plane to sensor (in um) 
scale =0;              % Scale for lateral distance between object and camera plane (in the wavevector direction)
N = 512;

% - Mismatches
mismatch_select = 2; % select on which parameter mismatches are added (0 for theta / 1 for phi / 2 fr zim)
if mismatch_select == 0
    dphi=0;dzim=0;
    dtheta=logspace(-1,1,10)/180*pi;
elseif mismatch_select == 1
    dtheta=0;dzim=0;
    dphi=logspace(-1,1,10)/180*pi;
elseif mismatch_select == 2
    dtheta=0; dphi=0;
    dzim=[0.01 0.025 0.05 0.075 0.1 0.15 0.2]*Lb;
end

% - Noise
sig_n = 2.5;                       % std of the gaussian noise will be nn*max(ynonoise(:)) 

% - Optim
MaxIter = 500;                          % Maximum number of iterations
cv_tol = 1e-8;                          % Tolerance for stopping criteria
RegParam =5e-8*sig_n.^2;                % Regularization (tuning) parameter 
RegType = 1;                            % Regularization type 
nonneg=1;
alt_iter = 1;                           % Number of alternate iterations
OTFatt=0;

%% Precomputations
% - Object
im = double(loadtiff('phantom.tif')); im=im/max(im(:));

% - Incident wave
theta = (theta_deg*pi)/180;        % Conversion in radian
phi = (phi_deg*pi)/180;            % Conversion in radian
k0 = 2*pi/Lb;                      % Wavenumber
ki = [sin(theta)*cos(phi),sin(theta)*sin(phi),-cos(theta)]*k0;    % Incident wavevector

% - Camera
lat_disp = (zim'./tan(theta))*scale;       % Lateral distance between object and camera plane (in the wavevector direction)
cam_shift=ki/norm(ki);
cam_shift=lat_disp*cam_shift/norm(cam_shift);  
px_sz=cam_side/N;

%% Data generation
% - Fwd Model
H = BuildFwd(cam_side, N, zim, ki, cam_shift, k0,0);
% - Noiseless data
ynonoise=H*im;
% - Noisy data
y=ynonoise+randn(size(ynonoise))*sig_n;

% - Displays
figure
subplot(2,2,1); imagesc(ynonoise); axis image; title('Noiseless data');
subplot(2,2,2); imagesc(log(1+abs(fftshift(fft2(ynonoise))))); axis image; title('FFT noiseless data');
subplot(2,2,3); imagesc(y); axis image; title('Noisy data');
subplot(2,2,4); imagesc(log(1+abs(fftshift(fft2(y))))); axis image; title('FFT noisy data');
colormap(flipud(gray))

%% Reconstruction
for ii=1:max([length(dtheta),length(dphi),length(dzim)])
    % -- Reconstruction Full model
    xFull = ReconsFullModel(cam_side,Lb,theta+dtheta(min(ii,length(dtheta))),phi+dphi(min(ii,length(dphi))),zim+dzim(min(ii,length(dzim))),lat_disp,RegType,RegParam,nonneg,MaxIter,alt_iter,cv_tol,y,0);
    [xJointDeconv,xc,xs] = JointDeconv(cam_side,Lb,theta+dtheta(min(ii,length(dtheta))),phi+dphi(min(ii,length(dphi))),zim+dzim(min(ii,length(dzim))),lat_disp,RegType,RegParam,MaxIter,alt_iter,cv_tol,OTFatt,y,0);

    % -- Compute metrics
    % SNR
    psnrFull(ii) = OptSNR(im,xFull);
    psnrJointDeconv(ii) = OptSNR(im,xJointDeconv);

    % FRC
    p.beta=1;
    freq=linspace(0,1/(px_sz*2),N/2+1);
    frc = FSC(im,xFull,p);
    idx = find(frc.frc<frc.T_hbit,1,'first');
    resFull(ii) = 1/freq(idx);
    frc = FSC(im,xJointDeconv,p);
    idx = find(frc.frc<frc.T_hbit,1,'first');
    resJointDeconv(ii) = 1/freq(idx);

    % -- Displays
    disp(['dtheta: ',num2str(dtheta(min(ii,length(dtheta)))/pi*180),' | dphi: ',num2str(dphi(min(ii,length(dphi)))/pi*180),' | dzim (in Lb): ',num2str(dzim(min(ii,length(dzim)))/Lb),' | PSNR Full: ',num2str(psnrFull(ii)),' | PSNR JointDeconv: ',num2str(psnrJointDeconv(ii)),' | Res Full: ',num2str(resFull(ii)),' | Res JointDeconv: ',num2str(resJointDeconv(ii))]);
    
    if max([length(dtheta),length(dphi),length(dzim)])>1
        figure(2);clf;hold all;grid;
        if length(dtheta)>1
            xx = dtheta/pi*180;
            xlabel('Error on $\theta$ ($^\circ$)','Interpreter','Latex');
        elseif length(dphi)>1
            xx = dphi/pi*180;
            xlabel('Error on $\phi$ ($^\circ$)','Interpreter','Latex');
        elseif length(dzim)>1
            xx = dzim/Lb;
            xlabel('Error on $z_{\mathrm{im}}$','Interpreter','Latex');
        end
        plot(xx(1:length(resFull)),resFull,'d','Linewidth',2,'DisplayName', 'FM');
        plot(xx(1:length(resJointDeconv)),resJointDeconv,'s','Linewidth',2,'DisplayName', 'JD');
        axis([xx(1) xx(end) 0.5*Lb Lb]);
        set(gca,'YTick',[0 0.25 0.5 0.75 1]*Lb);
        hAxes = get(gca,'YAxis');
        hAxes.TickLabelInterpreter = 'latex';
        if length(dzim)>1
            ll = {};
            for jj = 1:1:length(dzim)
                ll{jj} = [num2str(xx(jj)),'$$\lambda$$'];
            end
            hAxes = get(gca,'XAxis');
            hAxes.TickLabelInterpreter = 'latex';
            set(gca,'XTick',xx(1:1:end));
            set(gca,'XTickLabel',ll);
        else
            set(gca,'xscale','log');
        end
        set(gca,'YTickLabel',{'0','$$\frac{\lambda}{4}$$','$$\frac{\lambda}{2}$$','$$\frac{3\lambda}{4}$$','$$\lambda$$'});
        ylabel('$$\hat{r}_\mathrm{lim}$$','Interpreter','Latex');
        set(gca,'Fontsize',14);
        legend show; legend('Interpreter','Latex','Location','northeast');drawnow;

        figure(3);clf;hold all;grid;
        if length(dtheta)>1
            xx = dtheta/pi*180;
            xlabel('Error on $\theta$ ($^\circ$)','Interpreter','Latex');
        elseif length(dphi)>1
            xx = dphi/pi*180;
            xlabel('Error on $\phi$ ($^\circ$)','Interpreter','Latex');
        elseif length(dzim)>1
            xx = dzim/Lb;
            xlabel('Error on $z_{\mathrm{im}}$','Interpreter','Latex');
        end
        plot(xx(1:length(resFull)),psnrFull,'d','Linewidth',2,'DisplayName', 'FM');
        plot(xx(1:length(resJointDeconv)),psnrJointDeconv,'s','Linewidth',2,'DisplayName', 'JD');
        set(gca,'XLim',[xx(1) xx(end)]);
        if length(dzim)>1
            ll = {};
            for jj = 1:1:length(dzim)
                ll{jj} = [num2str(xx(jj)),'$$\lambda$$'];
            end
            hAxes = get(gca,'XAxis');
            hAxes.TickLabelInterpreter = 'latex';
            set(gca,'XTick',xx(1:1:end));
            set(gca,'XTickLabel',ll);
        else
            set(gca,'xscale','log');
        end
        hAxes = get(gca,'YAxis');
        hAxes.TickLabelInterpreter = 'latex';
        ylabel('PSNR','Interpreter','Latex');
        set(gca,'Fontsize',14);
        legend show; legend('Interpreter','Latex','Location','northeast');drawnow;
    end

    % -- Save
    imwrite(uint8((xFull-min(xFull(:)))/(max(xFull(:))-min(xFull(:)))*255),['xFull_dtheta',num2str(dtheta(min(ii,length(dtheta)))),'_dphi',num2str(dphi(min(ii,length(dphi)))),'_dzim',num2str(dzim(min(ii,length(dzim)))),'.png']);
    imwrite(uint8((xJointDeconv-min(xJointDeconv(:)))/(max(xJointDeconv(:))-min(xJointDeconv(:)))*255),['xJointDeconv_dtheta',num2str(dtheta(min(ii,length(dtheta)))),'_dphi',num2str(dphi(min(ii,length(dphi)))),'_dzim',num2str(dzim(min(ii,length(dzim)))),'.png']);
end


figure
subplot(2,2,1); imagesc(xFull); axis image; title('Full Model');
subplot(2,2,3); imagesc(log(1+abs(fftshift(fft2(xFull))))); axis image; title('FFT Full Model');
subplot(2,2,2); imagesc(xJointDeconv); axis image; title('Joint Deconv');
subplot(2,2,4); imagesc(log(1+abs(fftshift(fft2(xJointDeconv))))); axis image; title('FFT Joint Deconv');
colormap(flipud(gray))