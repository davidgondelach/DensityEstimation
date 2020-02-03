function [xf_mee] = propagateState_MeeBcRom(x0_mee,t0,tf,AC,BC,SWinputs,r,nop,svs,F_U,M_U,maxAtmAlt,et0,jdate0)
%PROPAGATESTATE_MEEBCROM - Propagate objects and ROM density
% Convert state in modified equinoctial elements to Cartesian coordinates
% propagate states and reduced-order density and convert Cartesian states
% back to modified equinoctial elements.
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

mu = 398600.4415;

xx_pv = x0_mee;
for k = 1:nop
    for j=1:size(x0_mee,2)
        [pos,vel] = ep2pv(x0_mee((k-1)*svs+1:(k-1)*svs+6,j),mu);
        xx_pv((k-1)*svs+1:(k-1)*svs+3,j) = pos;
        xx_pv((k-1)*svs+4:(k-1)*svs+6,j) = vel;
    end
end

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events', @(t,x) isdecayed(t,x,nop*size(x0_mee,2),svs));
[~,xf_out]=ode113(@(t,x) computeDerivative_PosVelBcRom(t,x,AC,BC,SWinputs,r,nop,svs,F_U,M_U,maxAtmAlt,et0,jdate0),[t0 tf],xx_pv,opts);
xf_pv = reshape(xf_out(end,:)',nop*svs+r,[]);

xf_mee = xf_pv;
for k = 1:nop
    for j=1:size(xf_pv,2)
        pos = xf_pv(svs*(k-1)+1:svs*(k-1)+3,j);
        vel = xf_pv(svs*(k-1)+4:svs*(k-1)+6,j);
        xf_mee((k-1)*svs+1:(k-1)*svs+6,j) = pv2ep(pos,vel,mu)';
    end
    % Make sure the difference in true longitude L of the sigma points wrt
    % the nominal true longitude L0 is minimal (i.e. L-L0 <= pi)
    % If the nominal true longitude L0 is close to pi (i.e. pi/2<L0 or L0<-pi/2) then wrap 
    % all L to [0,2pi] domain, so all difference in L <=pi (by default L is
    % on [-pi,pi] domain).
    if xf_mee((k-1)*svs+6,1) > pi/2 || xf_mee((k-1)*svs+6,1) < -pi/2
        xf_mee((k-1)*svs+6,:) = wrapTo2Pi(xf_mee((k-1)*svs+6,:));
    end
end

end

%------------- END OF CODE --------------
