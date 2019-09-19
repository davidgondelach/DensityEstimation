% Calculation of transformation matrix from Cartesian to RTN coordinates
function [cart2rtnMatrix] = computeCart2RTNMatrix(rrCart, vvCart)

% input: Cartesian position vector
%        Cartesain velocity vector
% output: Transformation matrix from Cartesian to RTN coordinates

iir = rrCart / norm(rrCart);
iih = cross( rrCart, vvCart );
iih = iih / norm(iih);
iit = cross( iih, iir );
iit = iit / norm(iit);
cart2rtnMatrix = [iir, iit, iih]';

end