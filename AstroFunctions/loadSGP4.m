function [ ] = loadSGP4()
%LOADSGP4 - Load SGP4 constants
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 03-Oct-2019

global tumin mu radiusearthkm xke j2 j3 j4 j3oj2 opsmode whichconst
opsmode = 'i';
whichconst = 72;
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );

end

