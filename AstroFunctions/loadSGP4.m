function [ ] = loadSGP4()
%LOADSGP4 - Load SGP4 constants
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

global tumin mu radiusearthkm xke j2 j3 j4 j3oj2 opsmode whichconst
opsmode = 'i';
whichconst = 72;
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );

end

