%--------------------------------------------------------------------Air/Ag
clear all;

c       = 3*10^8;
nm      = 10^(-9);
lambda0 = 575*nm;
b       = 0.1*lambda0;
a       = 0.93*b;
f       = c/lambda0;
e1          = 1;
[e2, flag2] = SilverComplexPermittivity(f);

xmin = -0.2001*lambda0;
xmax = +0.2*lambda0;
ymin = -0.2*lambda0;
ymax = +0.2001*lambda0;

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
grid;

%--------------------------------------------------------------------Air/Al
clear all;

c       = 3*10^8;
nm      = 10^(-9);
lambda0 = 495*nm;
b       = 0.167*lambda0;
a       = 0.95*b;
f       = c/lambda0;
e1          = 1;
[e2, flag2] = AluminiumComplexPermittivity(f);

xmin = -0.2001*lambda0;
xmax = +0.2*lambda0;
ymin = -0.2*lambda0;
ymax = +0.2001*lambda0;

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
grid;

%--------------------------------------------------------------------Air/Au
clear all;

c       = 3*10^8;
nm      = 10^(-9);
lambda0 = 515*nm;
b       = 0.08*lambda0;
a       = 0.61*b;
f       = c/lambda0;
e1          = 1;
[e2, flag2] = GoldComplexPermittivity(f);

xmin = -0.2001*lambda0;
xmax = +0.2*lambda0;
ymin = -0.2*lambda0;
ymax = +0.2001*lambda0;

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
grid;

%--------------------------------------------------------------------Air/Ti
clear all;

c       = 3*10^8;
nm      = 10^(-9);
lambda0 = 500*nm;
b       = 0.141*lambda0;
a       = 0.89*b;
f       = c/lambda0;
e1          = 1;
[e2, flag2] = TitaniumComplexPermittivity(f);

xmin = -0.2001*lambda0;
xmax = +0.2*lambda0;
ymin = -0.2*lambda0;
ymax = +0.2001*lambda0;

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
grid;