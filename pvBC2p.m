function [p] = pvBC2p(pvBC)
% Discard velocity and BC
p = pvBC(1:3,:);
end

