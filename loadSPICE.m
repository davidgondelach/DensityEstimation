function [ ] = loadSPICE(kernelpath)

% load standard kernels and reference frames
% Clear cspice memory
cspice_kclear;
% Load SPK, PCK, LSK kernels
cspice_furnsh( kernelpath );

end

