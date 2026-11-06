%--------------------------------------------------------------------------
% Script that reproduce the numerical evaluation of the resolution limit
% presented in Figures 5 and 6 of [1].
% 
% [1] Fourier Analysis of Interference Scanning Optical Probe Microscopy.
%     IEEE Transactions on Computational Imaging, vol. 11, pp. 1206-1216, 2025.
%     Emmanuel Soubies and Wolfgang Bacsa. 
%
% 2025 - Emmanuel Soubies (emmanuel.soubies@cnrs.fr)
%--------------------------------------------------------------------------
clear;close all; rng(1);
run('../Install')

%% Parameters
% - Incident wave
Lb = 0.5;                      % Wavelength (in um)
theta_deg = 60;                % Angle of the incident wave with z axis in degree
phi_deg = 0;                   % Angle of the incident wave with x axis in degree

% - Camera
cam_side = 90*Lb;              % Length of camera side (in um)
zim = [0.25 0.5 1 2 3]*Lb;     % Distance from image plane to sensor (in um) [If array -> loop]
scale =0;                      % Scale for lateral distance between object and camera plane (in the wavevector direction)
N = 512;

% - Noise
sig_n = [2.5 0.1];                      % std of the gaussian noise will be nn*max(ynonoise(:)) [If array -> loop]

% - Optim
MaxIter = 1000;                         % Maximum number of iterations
cv_tol = 1e-8;                          % Tolerance for stopping criteria
RegParam =5e-8*sig_n.^2;                % Regularization (tuning) parameter [If sig_n is an array, need to be of the same length]
RegType = 1;                            % Regularization type 
nonneg=0;
alt_iter = 1;                           % Number of alternate iterations
OTFatt=0;

%% Precomputations
isloop = length(sig_n)>1 || length(zim)>1;

% - Object
im = double(loadtiff('phantom.tif')); im=im/max(im(:));

% - Crops
if isloop
    crop1_i = 155:204;
    crop1_j = 220:269;
    crop2_i = 226:269;
    crop2_j = 352:395;
    figcrop1=figure;
    %subplot(length(sig_n)+1,length(zim),1); imagesc(im(crop1_i,crop1_j));axis image; axis off; colormap gray;title('GT');
    figcrop2=figure;
    %subplot(length(sig_n)+1,length(zim),1); imagesc(im(crop2_i,crop2_j));axis image; axis off; colormap gray;title('GT');
end

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

%% Loop
for ii = 1:length(sig_n)
    for jj = 1:length(zim)
        % -- Data generation
        % Fwd Model
        H = BuildFwd(cam_side, N, zim(jj), ki, cam_shift(jj,:), k0, 0);
        % Noiseless data
        ynonoise=H*im;
        % Noisy data
        y=ynonoise+randn(size(ynonoise))*sig_n(ii);

        if ~isloop
            figure
            subplot(2,2,1); imagesc(ynonoise); axis image; title('Noiseless data');
            subplot(2,2,2); imagesc(log(1+abs(fftshift(fft2(ynonoise))))); axis image; title('FFT noiseless data');
            subplot(2,2,3); imagesc(y); axis image; title('Noisy data');
            subplot(2,2,4); imagesc(log(1+abs(fftshift(fft2(y))))); axis image; title('FFT noisy data');
            colormap(flipud(gray))
        else
            figure(figcrop1)
            subplot(2*length(sig_n),length(zim),(ii-1)*2*(length(zim))+jj); imagesc(y(crop1_i,crop1_j));axis image; set(gca,'XTick',{});   set(gca,'YTick',{});colormap gray;
            if jj==1, ylabel([num2str(sig_n(ii)),'% noise']); end
            if ii==1, title(['$$z_{im} = $$',num2str(zim(jj)/Lb),'$$\lambda$$'],'Interpreter','latex'); end; drawnow;
            tmp = y(crop1_i,crop1_j);
            imwrite(uint8((tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)))*255),['y_Zoom1_sign_idx',num2str(ii),'_zim_idx',num2str(jj),'.png']);
            figure(figcrop2)
            subplot(2*length(sig_n),length(zim),(ii-1)*2*(length(zim))+jj); imagesc(y(crop2_i,crop2_j));axis image; set(gca,'XTick',{});   set(gca,'YTick',{}); colormap gray;
            if jj==1, ylabel([num2str(sig_n(ii)),'% noise']); end
            if ii==1, title(['$$z_{im} = $$',num2str(zim(jj)/Lb),'$$\lambda$$'],'Interpreter','latex'); end; drawnow;
            tmp = y(crop2_i,crop2_j);
            imwrite(uint8((tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)))*255),['y_Zoom2_sign_idx',num2str(ii),'_zim_idx',num2str(jj),'.png']);
        end

        % -- Reconstruction
        xest = ReconsFullModel(cam_side,Lb,theta,phi,zim(jj),lat_disp(jj),RegType,RegParam(ii),nonneg,MaxIter,alt_iter,cv_tol,y,0);
        figure(200); plot(mean(xest(63:134,227:260),1)); hold all;plot(mean(im(63:134,227:260),1))

        % -- Compute metrics
        % SNR
        snrTab(ii,jj) = OptSNR(im,xest);

        % FRC
        p.beta=1;
        freq=linspace(0,1/(px_sz*2),N/2+1);
        frc = FSC(im,xest,p);
        idx = find(frc.frc<frc.T_hbit,1,'first');
        res(ii,jj) = 1/freq(idx);

        % -- Plots
        if ~isloop
            figure;
            hold on
            plot(freq, frc.frc,'Linewidth',1.5, 'DisplayName', 'FRC')
            plot(freq, frc.T_hbit, 'Linewidth',1.5, 'DisplayName', '1/2 bit Threshold')
            plot([freq(idx) freq(idx)],[0,1],'--k', 'Linewidth',1.5, 'DisplayName', '$$\hat{\nu}_{\mathrm{lim}}$$');
            axis([0 1/(px_sz*2) 0 1]);
            legend show; grid
            legend('Interpreter','Latex');            
            set(gca,'XTick',[0 1/Lb 2/Lb 1/(px_sz*2)]);
            hAxes = get(gca,'XAxis');
            hAxes.TickLabelInterpreter = 'latex';
            hAxes = get(gca,'YAxis');
            hAxes.TickLabelInterpreter = 'latex';
            set(gca,'XTickLabel',{'0','$$k_0$$','$$2k_0$$','$$\frac{\nu_{\mathrm{s}}}{2}$$'});
            title('FRC Curve');
            set(gca,'Fontsize',14);

            figure;
            subplot(1,2,1); imagesc(xest); axis image; title('Full model');
            subplot(1,2,2); imagesc(log(1+abs(fftshift(fft2(xest))))); axis image; title('FFT full model');
            colormap(flipud(gray))           
        else
            figure(figcrop1)
            subplot(2*length(sig_n),length(zim),((ii-1)*2+1)*(length(zim))+jj); imagesc(xest(crop1_i,crop1_j));axis image; set(gca,'XTick',{});   set(gca,'YTick',{});colormap gray; drawnow;
            tmp = xest(crop1_i,crop1_j);
            imwrite(uint8((tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)))*255),['rec_Zoom1_sign_idx',num2str(ii),'_zim_idx',num2str(jj),'.png']);
            figure(figcrop2)
            subplot(2*length(sig_n),length(zim),((ii-1)*2+1)*(length(zim))+jj); imagesc(xest(crop2_i,crop2_j));axis image; set(gca,'XTick',{});   set(gca,'YTick',{});colormap gray; drawnow;
            tmp = xest(crop2_i,crop2_j);
            imwrite(uint8((tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)))*255),['rec_Zoom2_sign_idx',num2str(ii),'_zim_idx',num2str(jj),'.png']);
        end

        PlotResDiag(4,zim,res,sig_n,Lb,k0)
        disp(['sig_n: ',num2str(sig_n(ii)),' | zim: ',num2str(zim(jj)),' | SNR: ',num2str(snrTab(ii,jj)),' | Res: ',num2str(res(ii,jj))]);
    end
end

function PlotResDiag(fig,z,r,ss,Lb,k0)
    figure(fig);clf;hold all;grid;
    zim_val = logspace(-1,0.5,50);
    neff = [0.001 0.003 0.008 0.02 0.05 0.1 0.2];
    v_eff = @(neff,zim) sqrt((Lambert_W(pi^2*zim / (2*neff)) ./(zim)).^2 + k0^2) - k0;
    res = @(neff,zim) Lb ./ (1 + v_eff(neff,zim) /k0);
    for ne = neff
        res_val = res(ne,zim_val*Lb);
        plot(zim_val*Lb,res_val,'--','color',[.7 .7 .7],'HandleVisibility','off','Linewidth',1);
    end
    axis([0.17 zim_val(end) 0.25 1]*Lb);
    set(gca,'YTick',[0 0.25 0.5 0.75 1]*Lb);
    hAxes = get(gca,'YAxis');
    hAxes.TickLabelInterpreter = 'latex';
    set(gca,'YTickLabel',{'0','$$\frac{\lambda}{4}$$','$$\frac{\lambda}{2}$$','$$\frac{3\lambda}{4}$$','$$\lambda$$'});
    set(gca,'XTick',[0.17 0.5 1 2 3]*Lb);
    hAxes = get(gca,'XAxis');
    hAxes.TickLabelInterpreter = 'latex';
    set(gca,'XTickLabel',{'$$\frac{\lambda}{6}$$','$$\frac{\lambda}{2}$$','$$\lambda$$','$$2 \lambda$$','$$3 \lambda$$'});
    xlabel('$$z_\mathrm{im}$$','Interpreter','Latex');ylabel('$$\hat{r}_\mathrm{lim}$$','Interpreter','Latex');
    set(gca,'Fontsize',14);
    dd={'d','s','*','^','o'};
    for ii=1:size(r,1)
        plot(z(1:size(r,2)),r(ii,:),dd{ii},'Linewidth',2,'DisplayName', ['$$\hat{r}_\mathrm{lim}$$ (',num2str(ss(ii)),'\% noise)']);
    end
    legend show; legend('Interpreter','Latex','Location','southeast');drawnow;
end
