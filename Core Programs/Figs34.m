%--------------------------------------------------------------------SILVER
clear all;

nm          = 10^(-9);
c           = 3*10^8;
lambda0     = 575*nm;
k0          = 2*pi/lambda0;
f           = c/lambda0;

e1 = 1;
[e2, flag2] = SilverComplexPermittivity(f);

b = 0.1*lambda0;
a = 0.93*b;

dmin = -0.3*lambda0;
dmax = +0.3*lambda0;

[Hzx1, xo1, Hzy1, yo1] =...
    TMFieldDistributions1D(dmin, dmax, a, b, e1, e2, lambda0);

figure;
hold;
plot(xo1/lambda0, real(Hzx1), '-b', 'Linewidth', 2);
plot(yo1/lambda0, real(Hzy1), '-r', 'Linewidth', 3);
plot([-b/lambda0 -b/lambda0],  [-3, 3], 'g--', 'Linewidth', 1);
plot([-a/lambda0 -a/lambda0],  [-3, 3], 'g--', 'Linewidth', 1);
plot(-[-b/lambda0 -b/lambda0], [-3, 3], 'g--', 'Linewidth', 1);
plot(-[-a/lambda0 -a/lambda0], [-3, 3], 'g--', 'Linewidth', 1);
set(gca, 'fontsize', 16, 'fontname', 'times');
xlabel('d/\lambda_0');
ylabel('Re[H_z]');
legend('d=x (y=0)', 'd=y (x=0)');
ylim([-0.5, 2]);
grid;

% figure;
% hold;
% plot(xo1/lambda0, abs(Hzx1), '-b', 'Linewidth', 2);
% plot(yo1/lambda0, abs(Hzy1), '-r', 'Linewidth', 3);
% plot([-b/lambda0 -b/lambda0],  [0, 3], 'g--', 'Linewidth', 1);
% plot([-a/lambda0 -a/lambda0],  [0, 3], 'g--', 'Linewidth', 1);
% plot(-[-b/lambda0 -b/lambda0], [0, 3], 'g--', 'Linewidth', 1);
% plot(-[-a/lambda0 -a/lambda0], [0, 3], 'g--', 'Linewidth', 1);
% set(gca, 'fontsize', 16, 'fontname', 'times');
% xlabel('d/\lambda_0');
% ylabel('|H_z|');
% legend('d=x (y=0)', 'd=y (x=0)');
% % ylim([-0.5, 2]);
% grid;

%------------------------------------------------------------------ALUMINUM
clear all;

nm          = 10^(-9);
c           = 3*10^8;
lambda0     = 495*nm;
k0          = 2*pi/lambda0;
f           = c/lambda0;

e1 = 1;
[e2, flag2] = AluminiumComplexPermittivity(f);

b = 0.167*lambda0;
a = 0.95*b;

dmin = -0.3*lambda0;
dmax = +0.3*lambda0;

[Hzx1, xo1, Hzy1, yo1] =...
    TMFieldDistributions1D(dmin, dmax, a, b, e1, e2, lambda0);

figure;
hold;
plot(xo1/lambda0, real(Hzx1), '-b', 'Linewidth', 2);
plot(yo1/lambda0, real(Hzy1), '-r', 'Linewidth', 3);
plot([-b/lambda0 -b/lambda0],  [-3, 3], 'g--', 'Linewidth', 1);
plot([-a/lambda0 -a/lambda0],  [-3, 3], 'g--', 'Linewidth', 1);
plot(-[-b/lambda0 -b/lambda0], [-3, 3], 'g--', 'Linewidth', 1);
plot(-[-a/lambda0 -a/lambda0], [-3, 3], 'g--', 'Linewidth', 1);
set(gca, 'fontsize', 16, 'fontname', 'times');
xlabel('d/\lambda_0');
ylabel('Re[H_z]');
legend('d=x (y=0)', 'd=y (x=0)');
ylim([-0.5, 2]);
grid;

% figure;
% hold;
% plot(xo1/lambda0, abs(Hzx1), '-b', 'Linewidth', 2);
% plot(yo1/lambda0, abs(Hzy1), '-r', 'Linewidth', 3);
% plot([-b/lambda0 -b/lambda0],  [0, 3], 'g--', 'Linewidth', 1);
% plot([-a/lambda0 -a/lambda0],  [0, 3], 'g--', 'Linewidth', 1);
% plot(-[-b/lambda0 -b/lambda0], [0, 3], 'g--', 'Linewidth', 1);
% plot(-[-a/lambda0 -a/lambda0], [0, 3], 'g--', 'Linewidth', 1);
% set(gca, 'fontsize', 16, 'fontname', 'times');
% xlabel('d/\lambda_0');
% ylabel('|H_z|');
% legend('d=x (y=0)', 'd=y (x=0)');
% % ylim([-0.5, 2]);
% grid;