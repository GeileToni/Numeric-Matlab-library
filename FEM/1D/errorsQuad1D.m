% author: Mauro Morini
% last modified: 28.10.23
function [L2err, H1err] = errorsQuad1D(T, x, uh, dudx, u)
% calculates the error in L^2 and H^1 norm for quadratic FE with
% uh: numerical solution
% u: exact solution
% dudx: exact 1-st derivative of u

r = size(T, 2);

L2err = 0;
H1err = 0;

% weights QF
w = [5/9 8/9 5/9];

% nodes of QF
y = [-sqrt(3/5), 0, sqrt(3/5)];

% shape functions
N = {@(xi) 1/2*(xi.^2 - xi); @(xi) 1-xi.^2; @(xi) 1/2*(xi.^2 + xi)};
dN = {@(xi)1/2*(2*xi-1), @(xi) -2*xi, @(xi) 1/2*(1+2*xi)};

% iterate over elements
for i = 1:size(T, 1)

    % element
    K = x(T(i, :));

    % element length
    h = abs(K(1) - K(end));

    % middle of element
    m = (K(end) + K(1))/2;

    % element map
    F = @(xi) m + xi*h/2;
    
    uhVal = 0;
    duhVal = 0;
    for p = 1:r
        uhVal = uhVal + uh(T(i, p))*N{p}(y);
        duhVal = duhVal + uh(T(i, p))*dN{p}(y);
    end

    L2err = L2err + h/2*real(u(F(y)) - uhVal).^2*w.';
    H1err = H1err + h/2*real(dudx(F(y)) - 2/h*duhVal).^2*w.';
  
end
H1err = sqrt(H1err + L2err);
L2err = sqrt(L2err);

end