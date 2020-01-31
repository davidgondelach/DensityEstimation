function [ ] = loadSPICE(kernelpath)
%LOADSPICE - Load SPICE kernels and ephemerides
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 03-Oct-2019


% load standard kernels and reference frames
% Clear cspice memory
cspice_kclear;
% Load SPK, PCK, LSK kernels
cspice_furnsh( kernelpath );

end

