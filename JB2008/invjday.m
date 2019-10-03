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
% invjday: converts Julian date to Gregorian (calendar) date
%
% Input:
%  jd 		Julian date
%
% outputs:
%  month	calendar month [1-12]
%  day		calendar day [1-31]
%  year		calendar year [yyyy]
%  hr		hour
%  min		minute
%  sec		second
%
%  note: day may include fractional part
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
function [year, month, day, hr, min, sec] = invjday(jd)

z = fix(jd + .5);
fday = jd + .5 - z;

if (fday < 0)
   fday = fday + 1;
   z = z - 1;
end

if (z < 2299161)
   a = z;
else
   alpha = floor((z - 1867216.25) / 36524.25);
   a = z + 1 + alpha - floor(alpha / 4);
end
 
b = a + 1524;
c = fix((b - 122.1) / 365.25);
d = fix(365.25 * c);
e = fix((b - d) / 30.6001);
day = b - d - fix(30.6001 * e) + fday;
 
if (e < 14)
   month = e - 1;
else
   month = e - 13;
end
 
if (month > 2)
   year = c - 4716;
else
   year = c - 4715;
end

hr = abs(day-floor(day))*24;
min = abs(hr-floor(hr))*60;
sec = abs(min-floor(min))*60;

day = floor(day);
hr = floor(hr);
min = floor(min);

