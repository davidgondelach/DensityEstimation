function [Inputs] = Comp_Inputs_TIEGCM(jd0,jdf,TIEGCMSWdata)
%Comp_Inputs_TIEGCM - Compute space weather inputs for ROM-TIEGCM model

% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Sep 2019; Last revision: 24-Sep-2019

%------------- BEGIN CODE --------------

tt = jd0:1/24:jdf;
nofPoints = length(tt);

jd19970101 = 2450449.5;

Inputs = zeros(6,nofPoints);
for i=1:nofPoints
    jdate = tt(i);
    row = round((jdate-jd19970101)*24);
    % [jdate; doy; hour; f107d; f107a; Kp]
    Inputs(1,i) = jdate;
    Inputs(2:6,i) = TIEGCMSWdata(row,:)';
end
% Add future (now+1hr) values F10 and Kp
Inputs(7,1:end-1) = Inputs(4,2:end); % F10 (now+1hr)
Inputs(8,1:end-1) = Inputs(6,2:end); % Kp (now+1hr)
Inputs(7,end) = TIEGCMSWdata(row+1,3)'; % F10 (now+1hr)
Inputs(8,end) = TIEGCMSWdata(row+1,5)'; % Kp (now+1hr)

% Add quadratic Kp
Inputs(9,:) = Inputs(6,:).^2; % Kp^2 (now)
Inputs(10,:) = Inputs(8,:).^2; % Kp^2 (now+1hr)
% Add mixed terms F10*Kp
Inputs(11,:) = Inputs(6,:).*Inputs(4,:); % Kp*F10 (now)
Inputs(12,:) = Inputs(8,:).*Inputs(7,:); % Kp*F10 (now+1hr)

end

%------------- END OF CODE --------------
