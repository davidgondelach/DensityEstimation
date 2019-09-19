function [xf_pv] = propagateStatePosVel_FullGravDrag_ROM(x0_pv,time,et0,romStateTime,r,F_U,M_U,varargin)

% noo = 1;
% svs = 7;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
% [~,xf_pv]=ode113(@(t,x) Propagation_FullGravDrag_NRLMSISE(t,x,et0,noo,svs),time,x0_pv,opts);
[~,xf_pv]=ode113(@(t,x) Propagation_FullGravDrag_ROM(t,x,et0,romStateTime,r,F_U,M_U,varargin),time,x0_pv,opts);
xf_pv = xf_pv';

end
