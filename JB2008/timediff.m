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
% timediff: Time differences [s]
% 
% Last modified:   2018/01/27   M. Mahooti
% 
%--------------------------------------------------------------------------
function [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC)

TT_TAI  = +32.184;          % TT-TAI time difference [s]

GPS_TAI = -19.0;            % GPS-TAI time difference [s]

TT_GPS  =  TT_TAI-GPS_TAI;  % TT-GPS time difference [s]

TAI_GPS = -GPS_TAI;         % TAI-GPS time difference [s]

UT1_TAI = UT1_UTC-TAI_UTC;  % UT1-TAI time difference [s]

UTC_TAI = -TAI_UTC;         % UTC-TAI time difference [s]
  
UTC_GPS = UTC_TAI-GPS_TAI;  % UTC_GPS time difference [s]

UT1_GPS = UT1_TAI-GPS_TAI;  % UT1-GPS time difference [s]

TT_UTC  = TT_TAI-UTC_TAI;   %  TT-UTC time difference [s]

GPS_UTC = GPS_TAI-UTC_TAI;  % GPS-UTC time difference [s]

