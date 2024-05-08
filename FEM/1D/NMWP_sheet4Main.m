% author: Mauro Morini
% last modified: 22.03.24
clc;clear;close all;

% Initalizations
H = 2*pi*[40, 80, 160, 320].^(-1);
k = 4;
dirichlet = false;
error = zeros(4, length(H));
i = 1;

% functions
c = @(x) ones(size(x));
u_exact = @(x) exp(1i*k*x);
gTest = @(x) ones(size(x));
g = @(x) (x == 0);
dudx = @(x) 1i*k*exp(1i*k*x);

for i = 1:length(H)
    % grid
    h = 1/5;
    x = 0:h:(2*pi);
    n = length(x);
    Tlin = [(1:(n-1))', (2:n)'];
    Tquad = [(1:(n-2))', (2:(n-1))', (3:(n))'];
    int = true(size(x));
    int(1) = false;
    
    % get matrices
    Alin = StiffnessMatrix1D(Tlin,x,c);
    Aquad = StiffnessMatrix1D(Tquad,x,c);
    Mlin = massMatrix1D(x, Tlin, c);
    Mquad = massMatrix1D(x, Tquad, c);
    Blin = sparse(n,n);
    Blin(end,end) = 1i*k;
    Bquad = sparse(n,n);
    Bquad(end,end) = 1i*k;
    Uexact = real(u_exact(x));
    Ulin = zeros(n, 1);
    Uquad = zeros(n, 1);
    
    % calculate test Operator for double dirichlet b.c. 
    Klin = (-Alin + k^2*Mlin);
    Kquad = (-Aquad + k^2*Mquad);
    
    % solve pure dirichlet system
    if dirichlet
        intTest = int;
        intTest(end) = false;
        GTest = gTest(x');
        RHSLinTest = -Klin*GTest;
        RHSQuadTest = -Kquad*GTest;
        Ulin(intTest) = Klin(intTest, intTest)\RHSLinTest(intTest);
        Uquad(intTest) = Kquad(intTest, intTest)\RHSQuadTest(intTest);
        Uquad = Uquad + GTest;
        Ulin = Ulin + GTest;
    
        % calculate error
        [L2err, H1err] = errorsLinear1D(Tlin, x, Ulin, dudx, u_exact);
        error(1, i) = L2err;
        error(2, i) = H1err;
        [L2err, H1err] = errorsQuad1D(Tquad, x, Uquad, dudx, u_exact);
        error(3, i) = L2err;
        error(4, i) = H1err;
        
        figure(1)
        plot(x, Uquad, x, Ulin, x, Uexact, 'x')
        legend("Uquad", "Ulin", "Uexact")
    end
    
    % solve system
    if ~dirichlet
        Kquad = Kquad + Bquad;
        Klin = Klin + Blin;
        G = g(x');
        LHSlin = -Klin*G;
        LHSquad = -Kquad*G;
        Ulin(int) = Klin(int, int)\LHSlin(int);
        Uquad(int) = Kquad(int, int)\LHSquad(int);
        Ulin = real(Ulin + G);
        Uquad = real(Uquad + G);
    
        % calculate error
        [L2err, H1err] = errorsLinear1D(Tlin, x, Ulin, dudx, u_exact);
        error(1, i) = L2err;
        error(2, i) = H1err;
        [L2err, H1err] = errorsQuad1D(Tquad, x, Uquad, dudx, u_exact);
        error(3, i) = L2err;
        error(4, i) = H1err;
        
        % plot
        figure(1)
        plot(x, Uquad, x, Ulin, x, Uexact, 'x')
        legend("Uquad", "Ulin", "Uexact")
    end
end
error = real(error);
figure(2)
loglog(H, error(1,:),H, error(2,:),H, error(3,:),H, error(4,:), H, H.^2, '--', H, H, '--')
legend("L2error linear", "H1error linear", "L2error quadratic", "H1error quadratic", "H^2", "H")
