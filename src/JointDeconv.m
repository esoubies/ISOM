function [x,xc,xs] = JointDeconv(cam_side,Lb,theta,phi,h,lat_disp,RegType,RegParam,MaxIter,alt_iter,cv_tol,OTFatt,y,displ) 
%--------------------------------------------------------------------------
% Joint deconvolution approach described in Section IV.B.2 of [1]
%
% Inputs : cam_side  -> Length of camera side (in nm)
%          Lb        -> Wavelength (in um)
%          theta     -> Angle of the incident wave with z axis in radian
%          phi       -> Angle of the incident wave with x axis in radian
%          h         -> Distance from image plane to sensor (in nm)
%          lat_disp  -> Lateral shift of the camera center 
%          RegType   -> regularization choice (1: order-0 Tikhonov, 2: order-1 Tikhonov, 3: smoothed-TV)
%          RegParam  -> regularization parameter
%          MaxIter   -> Max number of iterations of inner solver (optim cos and sin parts)
%          alt_iter  -> Number of alternate iter (between optim offset and cos/sin parts optim)
%          cv_tol    -> stopping criteria threshold (relative cost and relative step)
%          OTFatt    -> boolean true if OTF attenuation is activated
%          y         -> Measurements image
%          disp      -> boolean true to display results
%
% Output : x         -> reconstructed image
%          xc        -> cos part of reconstructed image
%          xs        -> sin part of reconstructed image
%
% [1] Fourier Analysis of Interference Scanning Optical Probe Microscopy.
%     IEEE Transactions on Computational Imaging, vol. 11, pp. 1206-1216, 2025.
%     Emmanuel Soubies and Wolfgang Bacsa. 
%
% 2025 - Emmanuel Soubies (emmanuel.soubies@cnrs.fr)
%--------------------------------------------------------------------------

if nargin < 14
    displ=1;
end

%% Pre-computations
% - Physical parameters
M = size(y,1); N=M;
k0 = 2*pi/Lb;                                                                     
ki = [sin(theta)*cos(phi),sin(theta)*sin(phi),-cos(theta)]*k0;                    % Incident wavevector
cam_shift=ki/norm(ki);
cam_shift=cam_shift/norm(cam_shift)*lat_disp;  

% - Data normalization
m = max(mean(y(:)),0);
y=y/norm(y(:)-m);

%% Definition of operators and costs
% - Enlarged grid (to account for zero-padding)
rj=linspace(-cam_side,cam_side-cam_side/N,2*N);
[X,Y]=meshgrid(rj);Y=Y(end:-1:1,:);

Xpsf=X+cam_shift(1);
Ypsf=Y+cam_shift(2);
norm_r=sqrt(Xpsf.^2+Ypsf.^2+h^2);
hhs = -sin((ki(1)*Xpsf + ki(2)*Ypsf) - k0*norm_r) ./ norm_r;
hhc = cos((ki(1)*Xpsf + ki(2)*Ypsf) - k0*norm_r) ./ norm_r;
nn = max(norm(hhs(:),1),norm(hhc(:),1));
hhs=hhs/nn;hhc=hhc/nn;

% - Operators
HH = LinOpConv(fft2(ifftshift(cat(3,hhc,hhs))),1,[1,2]);
indxmin = [1,1,1]; indxmax = [N,N,2];
P = LinOpSelectorPatch(HH.sizein,indxmin,indxmax); 
SS = LinOpSum(P.sizeout,3);
H = SS*P*HH*P';  
Sq = OpEWSquaredMagnitude(P.sizeout);
Sqr = OpEWSqrt(SS.sizeout);
% - Regularizer
if RegType == 1
    R=RegParam*CostL2(SS.sizeout)*Sqr*SS*Sq;
elseif RegType == 2
    G = LinOpGrad(SS.sizeout);
    R=RegParam*CostL2(G.sizeout)*G*Sqr*SS*Sq;
elseif RegType == 3
    %G = LinOpGrad(SS.sizeout);
    R=RegParam*CostHyperBolic(Sqr.sizeout,1e-5,[])*Sqr*SS*Sq;
end

% Generate OTF att
if OTFatt
    lv = ((1:M) - floor(M/2)-1)/M;
    lh = ((1:M) - floor(M/2)-1)/M;
    [X,Y]=meshgrid(lh,lv);
    [~,rho]=cart2pol(X,Y);

    width=0.05;
    str=0.99;
    OTFatt=ifftshift((1-str*exp(-rho.^2/(2*(width)^2))));
    Hatt=LinOpConv('MTF', OTFatt,1,[1,2]);
end

%% Initialize optimization variables
m=max(mean(y(:)),0);
x = H'*(y-m);

%% Loop over outter iterations
if displ, fig=figure; end
for it=1:alt_iter
    % - Data Term
    if OTFatt
        D=CostL2([],y - m,Hatt)*H;
    else
        D=CostL2([],y - m)*H;
    end

    % - Optimization cos and sin parts
    VMLMB = OptiVMLMB(D+R);                                                       % Declare Forward Backward Splitting optimizer                                   % Descent gradient
    VMLMB.maxiter = MaxIter;
    VMLMB.OutOp = OutputOpti(1,round(MaxIter/10));
    VMLMB.ItUpOut = round(MaxIter/10);
    VMLMB.CvOp = TestCvgCombine('CostRelative',cv_tol,'SetRelative',cv_tol);
    VMLMB.run(x);

    % - Update x
    x=VMLMB.xopt;

    % - Update constant m
    diff = y-H*x;
    m = max(mean(diff(:)),0);

    % - Displays
    if displ
        figure(fig)
        yy = P*HH*P'*x;
        subplot(2,4,1);imagesc(x(:,:,1)); axis image; title('Cos part');
        subplot(2,4,2);imagesc(x(:,:,2)); axis image; title('Sin part');
        subplot(2,4,3);imagesc(sqrt(sum(x.^2,3))); axis image; title('Recons');
        subplot(2,4,5);imagesc(yy(:,:,1)); axis image; title('Cos part (y)');
        subplot(2,4,6);imagesc(yy(:,:,2)); axis image; title('Sin part (y)');
        subplot(2,4,7);imagesc(SS*yy+m); axis image; title('Recons (y)');
        subplot(2,4,8);imagesc(y); axis image; title('Data');colormap gray;
    end
end
xc = VMLMB.xopt(:,:,1);
xs = VMLMB.xopt(:,:,2);
x=sqrt(sum(VMLMB.xopt.^2,3));
end