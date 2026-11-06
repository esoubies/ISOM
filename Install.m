function Install()
%--------------------------------------------------------------------------
% Function Install()
%
% Ensures that all the necessray functions of ISOM and GlobalBioIm are
% in the MATLab path.
%
% 2025 - Emmanuel Soubies (emmanuel.soubies@cnrs.fr)
%--------------------------------------------------------------------------
addpath(pwd);
addpath(genpath([pwd,'/src']));     % Add source files to path
addpath(genpath([pwd,'/Utils'])); 

if exist('GlobalBioIm', 'dir')             % Check existance for GlobalBioIm in path
    GBIPath = getfield(what('GlobalBioIm'),'path');
    addpath(genpath(GBIPath))
else
    disp('GlobalBioIm not found in computer. Dowloading in current directory...');
    websave([pwd,'/GlobalBioIm'], 'https://github.com/Biomedical-Imaging-Group/GlobalBioIm/archive/refs/heads/master.zip');
    unzip('GlobalBioIm.zip')
    movefile('GlobalBioIm-master', 'GlobalBioIm')
    delete('GlobalBioIm.zip')
    run('GlobalBioIm/setGlobalBioImPath.m')
    savepath
end
