function [xf_pv] = propagateStatePosVel_FullGravDrag_NRLMSISE(x0_pv,time,et0,varargin)

noo = 1;
svs = 7;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
if nargin > 3
    [~,xf_pv]=ode113(@(t,x) Propagation_FullGravDrag_NRLMSISE(t,x,et0,noo,svs,varargin{1}),time,x0_pv,opts);
else
    [~,xf_pv]=ode113(@(t,x) Propagation_FullGravDrag_NRLMSISE(t,x,et0,noo,svs),time,x0_pv,opts);
end
xf_pv = xf_pv';

end
