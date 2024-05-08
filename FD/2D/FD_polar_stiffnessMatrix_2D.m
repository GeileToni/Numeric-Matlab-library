% author: Mauro Morini
% last modified: 10.04.24
function LHS = FD_polar_stiffnessMatrix_2D(ntheta, nr,r0,r1)
% calculates stiffness matrix for FD in 2D on a polar grid. i.e the FD
% discretization of d_rr u + 1/r*d_r u + 1/r^2*d_thetatheta u
% where the grid starts on the inner cirlce and ends ar the outer circle
%
% Inputs :
% dtheta: scalar number of grid points in theta direction (angle)
% dr : scalar number of grid points in r direction (length), does not
%       include points on inner circle
% r0 : scalar radius of inner obstacle
% r1 : scalar radius of artificial DtN boundary
%
% Output :
% LHS : (ntheta*nt,ntheta*nr) sparse matrix of the FD-stencil for the
%       laplacian in polar form

% Initializations
dtheta = 2*pi/ntheta;
dr = (r1 - r0)/nr;
ntot = ntheta*nr;
e = ones(ntheta, 1);

triDiagA = [e*(r1*dtheta)^(-2), e*(-2/dr^2 - 2/(r1*dtheta)^2), e*(r1*dtheta)^(-2)];
A = spdiags(triDiagA, -1:1, ntheta, ntheta);
A(1,ntheta) = 1/(r1*dtheta)^2;
A(ntheta,1) = 1/(r1*dtheta)^2;

%   diagonal block matrix on sub/super diagonal of LHS
Asub = spdiags(e*(1/dr^2 - (2*dr*r1)^(-1)),0,ntheta,ntheta);
Asuper = spdiags(e*(1/dr^2 + (2*dr*r1)^(-1)),0,ntheta,ntheta);

% fill LHS first row
LHS = sparse(ntot,ntot);
idxv = 1:ntheta;
LHS(idxv,idxv) = A;
LHS(idxv, ntheta + idxv) = Asuper;

% fill LHS second to second last row
for l = 1:(nr-2)
    LHS(l*ntheta + idxv,l*ntheta + idxv) = A;
    LHS(l*ntheta + idxv,(l-1)*ntheta + idxv) = Asub;
    LHS(l*ntheta + idxv,(l+1)*ntheta + idxv) = Asuper;
end

% fill LHS last row
LHS((nr-1)*ntheta + idxv, (nr-2)*ntheta + idxv) = Asub;
LHS((nr-1)*ntheta + idxv, (nr-1)*ntheta + idxv) = A;

end