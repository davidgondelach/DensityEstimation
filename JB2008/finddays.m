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

%---------------------------------------------------------------------------
%
%  finddays.m
%
%  this function finds the fractional days through a year given the year,
%  month, day, hour, minute and second.
%
%  Inputs:
%    year        - year                           1900 .. 2100
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - hour                           0 .. 23
%    min         - minute                         0 .. 59
%    sec         - second                         0.0 .. 59.999
%
%  Output:
%    days        - day of year plus fraction of a
%                  day                            days
%
%---------------------------------------------------------------------------
function days = finddays(year, month, day, hr, min, sec)

for i= 1:12
    lmonth(i) = 31;
    if i == 2
        lmonth(i)= 28;
    end
    if i == 4 || i == 6 || i == 9 || i == 11
        lmonth(i)= 30;
    end
end

if (rem (year,4) == 0)
    lmonth(2)= 29;
    if (rem (year,100) == 0) && (rem (year,400) ~= 0)
        lmonth(2)= 28;
    end
end

i = 1;
days = 0;
while (i < month) && ( i < 12 )
    days= days + lmonth(i);
    i= i + 1;
end

days= days + day + hr/24 + min/1440 + sec/86400;

