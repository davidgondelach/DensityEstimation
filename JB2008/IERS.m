% This file was downloaded on 11 June 2019 from:
% https://www.mathworks.com/matlabcentral/fileexchange/56163-jacchia-bowman-atmospheric-density-model
% 
% Copyright (c) 2018, Meysam Mahooti
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%--------------------------------------------------------------------------
%
% IERS: Management of IERS time and polar motion data
%  
% Last modified:   2018/02/01   M. Mahooti
% 
%--------------------------------------------------------------------------
function [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eop,Mjd_UTC,interp)

global const

if (nargin == 2)
   interp = 'n';
end

if (interp =='l')
    % linear interpolation
    mjd = (floor(Mjd_UTC));
    i = find(mjd==eop(4,:),1,'first');
    preeop = eop(:,i);
    nexteop = eop(:,i+1);
    mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
    fixf = mfme/1440;
    % Setting of IERS Earth rotation parameters
    % (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
    x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
    y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
	UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
    LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
    dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
    deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
    dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
    dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
    TAI_UTC = preeop(13);
	
    x_pole  = x_pole/const.Arcs;  % Pole coordinate [rad]
    y_pole  = y_pole/const.Arcs;  % Pole coordinate [rad]
    dpsi    = dpsi/const.Arcs;
    deps    = deps/const.Arcs;
    dx_pole = dx_pole/const.Arcs; % Pole coordinate [rad]
    dy_pole = dy_pole/const.Arcs; % Pole coordinate [rad]
elseif (interp =='n')    
    mjd = (floor(Mjd_UTC));
    i = find(mjd==eop(4,:),1,'first');
    eop = eop(:,i);
    % Setting of IERS Earth rotation parameters
    % (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
    x_pole  = eop(5)/const.Arcs;  % Pole coordinate [rad]
    y_pole  = eop(6)/const.Arcs;  % Pole coordinate [rad]
	UT1_UTC = eop(7);             % UT1-UTC time difference [s]
    LOD     = eop(8);             % Length of day [s]
    dpsi    = eop(9)/const.Arcs;
    deps    = eop(10)/const.Arcs;
    dx_pole = eop(11)/const.Arcs; % Pole coordinate [rad]
    dy_pole = eop(12)/const.Arcs; % Pole coordinate [rad]
	TAI_UTC = eop(13);            % TAI-UTC time difference [s]
end

