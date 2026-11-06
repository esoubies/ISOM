function xest = ReconsFullModel(cam_side,Lb,theta,phi,h,lat_disp,RegType,RegParam,nonneg,MaxIter,alt_iter,cv_tol,y,displ) 
%--------------------------------------------------------------------------
% Full Model reconstruction approach described in Section IV.B.1 of [1]
%
% Inputs : cam_side  -> Length of camera side (in nm)
%          Lb        -> Wavelength (in um)
%          theta     -> Angle of the incident wave with z axis in radian
%          phi       -> Angle of the incident wave with x axis in radian
%          h         -> Distance from image plane to sensor (in nm)
%          lat_disp  -> Lateral shift of the camera center 
%          RegType   -> regularization choice (1: order-0 Tikhonov, 2: order-1 Tikhonov, 3: smoothed-TV)
%          RegParam  -> regularization parameter
%          nonneg    -> boolean true is non-negativity constraint is activated
%          MaxIter   -> Max number of iterations of inner solver on the reconstructed image
%          alt_iter  -> Number of alternate iter (between optim offset and optim image)
%          cv_tol    -> stopping criteria threshold (relative cost and relative step)
%          y         -> Measurements image
%          disp      -> boolean true to display results
%
% Output : xest      -> reconstructed image
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
M = size(y,1);
k0 = 2*pi/Lb;                                                                     
ki = [sin(theta)*cos(phi),sin(theta)*sin(phi),-cos(theta)]*k0;                    % Incident wavevector
cam_shift=ki/norm(ki);
cam_shift=cam_shift/norm(cam_shift)*lat_disp;  

% - Data normalization
m = max(mean(y(:)),0);
y=y/norm(y(:)-m);

% - Build the fwd operator
N = M;
H = BuildFwd(cam_side, N, h, ki, cam_shift, k0,1);   

% - Initializations
m=max(mean(y(:)),0);
xest=H'*(y-m);

% Generate a normalized OTF
lv = ((1:M) - floor(M/2)-1)/M;
lh = ((1:M) - floor(M/2)-1)/M;
[X,Y]=meshgrid(lh,lv);
[~,rho]=cart2pol(X,Y);

% - Reconstruction
if displ
fig=figure;
end
for it=1:alt_iter   
    if RegType == 1
        R = RegParam*CostL2(H.sizein);
    elseif RegType == 2
        G = LinOpGrad(H.sizein);
        R = RegParam*CostL2(G.sizeout)*G;
    elseif RegType == 3
        G = LinOpGrad(H.sizein);
        R = RegParam*CostHyperBolic(G.sizeout,1e-5,3)*G;
    end
    L2=CostL2([],y - m);
    F = L2*H;
    if nonneg
        VMLMB = OptiVMLMB(F + R,0.,[]);
    else
        VMLMB = OptiVMLMB(F + R,[],[]);
    end
    VMLMB.maxiter = MaxIter;
    VMLMB.OutOp = OutputOpti(1,round(MaxIter/10));
    VMLMB.ItUpOut = round(MaxIter/10);
    VMLMB.CvOp = TestCvgCombine('CostRelative',cv_tol,'SetRelative',cv_tol);
    VMLMB.run(xest);                                                           % Minimization over x
    xest = VMLMB.xopt();
 
    diff = y-H*xest;
    m = max(mean(diff(:)),0);

    % - Displays
    if displ
    figure(fig)
    subplot(1,3,1);imagesc(y); axis image; title('Data');
    subplot(1,3,2);imagesc(xest); axis image; title('Recons');
    subplot(1,3,3);imagesc(H*xest); axis image; title('Recons (y)');colormap gray;
    end
end
end