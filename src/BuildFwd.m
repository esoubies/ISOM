function H = BuildFwd(cam_side, N, h, ki, cam_shift, k0,norm_conv) 
%--------------------------------------------------------------------------
% function H = BuildFwd(cam_side, N, h, ki, cam_shift, k0) 
%
% Generate the convolution operator H using convolution and diagonal
% operators, padding to mitigate edge effects and downsampling to equalize
% the dimensions the object plane and image plane. 
%
% Inputs : cam_side  -> Length of camera side (in nm)
%          N         -> Object grid of size N x N
%          h         -> Distance from image plane to sensor (in nm)
%          ki        -> Incident wavevector
%          cam_shift -> Lateral shift of the camera center 
%          k0        -> Wave number
%          norm_conv -> boolean true to normalize convolution kernels
%
% Output : H        -> Build operator for forward model
%
% 2025 - Emmanuel Soubies (emmanuel.soubies@cnrs.fr)
%--------------------------------------------------------------------------

%% Generate convolution operator
% - Enlarged grid (to account for zero-padding)
rj=linspace(-cam_side,cam_side-cam_side/N,2*N);
[X,Y]=meshgrid(rj);Y=Y(end:-1:1,:);

% - Convolution kernels
Xpsf=X+cam_shift(1);
Ypsf=Y+cam_shift(2);
norm_r=sqrt(Xpsf.^2+Ypsf.^2+h^2);
hf = cos((ki(1)*Xpsf + ki(2)*Ypsf) - k0*norm_r) ./ norm_r;
hg = sin((ki(1)*Xpsf + ki(2)*Ypsf) - k0*norm_r) ./ norm_r;
if norm_conv
    nn = max(norm(hf(:),1),norm(hg(:),1));
    hf=hf/nn;hg=hg/nn;
end

%hf=padarray(padarray(hf,ceil([N/2,N/2]),0,'pre'),floor([N/2,N/2]),0,'post');
%hg=padarray(padarray(hg,ceil([N/2,N/2]),0,'pre'),floor([N/2,N/2]),0,'post');

% - Convolution operators
Hf = LinOpConv(fft2(ifftshift(hf)));
Hg = LinOpConv(fft2(ifftshift(hg)));

%% Generate diagonal operators
% - Not enlarged grid
rj=linspace(-cam_side/2,cam_side/2-cam_side/N,N);
[X,Y]=meshgrid(rj);Y=Y(end:-1:1,:);

% - Diagonal values
df = cos(ki(1)*X+ki(2)*Y + ki(3)*h);
dg = sin(ki(1)*X+ki(2)*Y + ki(3)*h);

% - Diagonal operators
Df = LinOpDiag([],df);
Dg = LinOpDiag([],dg);

%% Padding and downsampling operators
sz = [2*N,2*N]; indxmin = [1,1]; indxmax = [N,N];
P = LinOpSelectorPatch(sz,indxmin,indxmax); 
H = P*(Hf*P'*Df-Hg*P'*Dg);                        %Complete forward model(to reconstruct on the same grid size)
end
