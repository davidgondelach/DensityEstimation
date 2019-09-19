function [xf_mee] = propagateStateMEE_FullGravDrag_NRLMSISE(x0_mee,t0,tf,et0)

mu = 398600.4415;
nop = 1;
svs = 7;

xx_pv = x0_mee;
for k = 1:nop
    for j=1:size(x0_mee,2)
        [pos,vel] = ep2pv(x0_mee((k-1)*svs+1:(k-1)*svs+6,j),mu);
        xx_pv((k-1)*svs+1:(k-1)*svs+3,j) = pos;
        xx_pv((k-1)*svs+4:(k-1)*svs+6,j) = vel;
    end
end

% opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'OutputFcn',@odewbar,'Events', @(t,x) isdecayed(t,x,nop*size(x0_mee,2),svs));
% [~,tempE]=ode113(@(t,x) Propagation11_ODE_Var_FullGrav_New(t,x,AC,BC,Inp2,r,nop,svs,F_U,M_U,et0,jdate0),[t0 tf],xx_pv,opts);
[tempE] = propagateStatePosVel_FullGravDrag_NRLMSISE(xx_pv,[t0,tf],et0);
xf_pv = reshape(tempE(:,end),nop*svs,[]);

xf_mee = xf_pv;
for k = 1:nop
    for j=1:size(xf_pv,2)
        pos = xf_pv(svs*(k-1)+1:svs*(k-1)+3,j);
        vel = xf_pv(svs*(k-1)+4:svs*(k-1)+6,j);
        xf_mee((k-1)*svs+1:(k-1)*svs+6,j) = pv2ep(pos,vel,mu)';
    end
end

end
