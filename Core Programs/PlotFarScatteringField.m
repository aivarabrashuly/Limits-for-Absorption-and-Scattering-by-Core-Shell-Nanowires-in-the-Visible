function [Fz1] = PlotFarScatteringField(C1)

PHI  = 1000;
phi1 = 2*pi*[0:PHI]/PHI;
N    = (length(C1)-1)/2;
n1   = [-N:N];

C2   = repmat(reshape(C1, [2*N+1, 1]),  [1, PHI+1]);
phi2 = repmat(reshape(phi1, [1, PHI+1]),[2*N+1, 1]);
n2   = repmat(reshape(n1, [2*N+1, 1]),  [1, PHI+1]);

Fz1 = sum(C2.*(1i.^n2).*exp(1i.*n2.*phi2));

figure;
polar(phi1, abs(Fz1).^2);
