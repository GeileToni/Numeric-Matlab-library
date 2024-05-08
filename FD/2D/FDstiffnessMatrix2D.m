% author: Mauro Morini  
% last modified: 26.04.24
function A = FDstiffnessMatrix2D(nx, ny)
% calculates not scaled FD stiffness matrix in 2D for laplacian with
% lexographical orientation (see Numerik der DGL)
% CAREFUL: don't forget to multiply 1/dx^2 to it
%
% Inputs: 
% nx: scalar number of interior points in x direction
% ny: scalar number of interior points in y direction
%
% Outputs:
% A: (nx*ny,nx*ny) non scaled stiffness matrix 

A = sparse(nx*ny, nx*ny);
Aloc = spdiags([-1, 4, -1], -1:1, nx,nx);
idxVec = 1:nx;
A(idxVec, idxVec) = Aloc;
A(idxVec, nx+idxVec) = speye(nx,nx)*(-1);
for j = 1:(ny-2)
    A(j*nx + idxVec,(j-1)*nx + idxVec) = speye(nx,nx)*(-1);
    A(j*nx + idxVec,(j+1)*nx + idxVec) = speye(nx,nx)*(-1);
    A(j*nx + idxVec,j*nx + idxVec) = Aloc;
end
A((ny-1)*nx + idxVec, (ny-1)*nx + idxVec) = Aloc;
A((ny-1)*nx + idxVec, (ny-2)*nx + idxVec) = speye(nx,nx)*(-1);

end