%--------------------------------------------------------------------Al/GaP
clear all;

c       = 3*10^8;
nm      = 10^(-9);
lambda0 = 400*nm;
b       = 0.089*lambda0;
a       = 0.32*b;
f       = c/lambda0;
[e1, flag1] = AluminiumComplexPermittivity(f);
[e2, flag2] = GalliumPhosphideComplexPermittivity(f);

xmin = -0.10001*lambda0;
xmax = +0.1*lambda0;
ymin = -0.1*lambda0;
ymax = +0.10001*lambda0;

[Ez2, x1, y1] =...
    TEFieldDistribution(xmin, xmax, ymin, ymax, a, b, e1, e2, lambda0);

PHI  = 1000;
phi1 = pi*[-PHI:PHI]/PHI;

figure;
contourf(x1/lambda0, y1/lambda0, abs(Ez2));
hold;
plot(a*cos(phi1)/lambda0, a*sin(phi1)/lambda0, '-w', 'LineWidth', 2.8);
plot(b*cos(phi1)/lambda0, b*sin(phi1)/lambda0, '-w', 'LineWidth', 2.8);
set(gca, 'fontsize', 16, 'fontname', 'times');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');
colorbar;
caxis([0, 2]);
grid;
% 
% k0 = 2*pi/lambda0;
% N  = 5;
% [CTE, CTM] = ScatteringCoefficients(e1, e2, a, b, k0, N);
% [PTEabs,  PTEscat,  PTMabs,  PTMscat,...
%           PTE0abs, PTE0scat, PTM0abs, PTM0scat] =...
%           NewScatteringAndAbsorbingPowers(CTE, CTM, k0, b);


%-------------------------------------------------------------------Ag/AlSb
clear all;

c       = 3*10^8;
nm      = 10^(-9);
lambda0 = 620*nm;
b       = 0.081*lambda0;
a       = 0.49*b;
f       = c/lambda0;
[e1, flag1] = SilverComplexPermittivity(f);
[e2, flag2] = AluminiumAntimonideComplexPermittivity(f);

xmin = -0.10001*lambda0;
xmax = +0.1*lambda0;
ymin = -0.1*lambda0;
ymax = +0.10001*lambda0;

[Hz2, x1, y1] =...
    TMFieldDistribution(xmin, xmax, ymin, ymax, a, b, e1, e2, lambda0);

PHI  = 1000;
phi1 = pi*[-PHI:PHI]/PHI;

figure;
contourf(x1/lambda0, y1/lambda0, abs(Hz2));
hold;
plot(a*cos(phi1)/lambda0, a*sin(phi1)/lambda0, '-w', 'LineWidth', 2.8);
plot(b*cos(phi1)/lambda0, b*sin(phi1)/lambda0, '-w', 'LineWidth', 2.8);
set(gca, 'fontsize', 16, 'fontname', 'times');
xlabel('x/\lambda_0');
ylabel('y/\lambda_0');
colorbar;
caxis([0, 14]);
grid;
% 
% k0 = 2*pi/lambda0;
% N  = 5;
% [CTE, CTM] = ScatteringCoefficients(e1, e2, a, b, k0, N);
% [PTEabs,  PTEscat,  PTMabs,  PTMscat,...
%           PTE0abs, PTE0scat, PTM0abs, PTM0scat] =...
%           NewScatteringAndAbsorbingPowers(CTE, CTM, k0, b);