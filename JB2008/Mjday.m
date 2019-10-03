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
% Mjday: Modified Julian Date from calendar date and time
%
% Inputs:
%   year, month, day, hour, min, sec
%
% output:
%   Mjd		Modified Julian Date
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function Mjd = Mjday(year, month, day, hour, min, sec)

if (nargin < 4)
    hour = 0;
    min  = 0;
    sec  = 0;
end

y = year;
m = month;
b = 0;
c = 0;

if (m <= 2)
   y = y - 1;
   m = m + 12;
end

if (y < 0)
   c = -.75;
end

% check for valid calendar date
if (year < 1582)
   % null
elseif (year > 1582)
   a = fix(y / 100);
   b = 2 - a + floor(a / 4);
elseif (month < 10)
   % null
elseif (month > 10)
   a = fix(y / 100);
   b = 2 - a + floor(a / 4);
elseif (day <= 4)
   % null
elseif (day > 14)
   a = fix(y / 100);
   b = 2 - a + floor(a / 4);
else
    fprintf('\n\n  this is an invalid calendar date!!\n');
    return
end

jd = fix(365.25 * y + c) + fix(30.6001 * (m + 1));
jd = jd + day + b + 1720994.5;
jd = jd + (hour+min/60+sec/3600)/24;
Mjd = jd - 2400000.5;

