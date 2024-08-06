clear all;

nm      = 10^(-9);
c       = 3*10^8;
N       = 6;

%-------------------------------------------------------------------SILVER
lambda0 = 575*nm;
b       = 0.1*lambda0;
a       = 0.93*b;

lambdamin = lambda0-200*nm;
lambdamax = lambda0+200*nm;
lambdapoints = 400;
lambda1   = lambdamin + (lambdamax-lambdamin)*[0:lambdapoints]/lambdapoints;

for lambdaindex = 0:lambdapoints
    lambda      = lambdamin + (lambdamax-lambdamin)*lambdaindex/lambdapoints;
    e1          = 1;
    e2          = SilverComplexPermittivity(c/lambda);
    [CTE, CTM]  = ScatteringCoefficients(e1, e2, a, b, 2*pi/lambda, N);
    [PTEabs1(1, lambdaindex+1),  PTEscat,...
     PTMabs1(1, lambdaindex+1),  PTMscat,...
    PTE0abs1(1, lambdaindex+1), PTE0scat,...
    PTM0abs1(1, lambdaindex+1), PTM0scat] =...
        NewScatteringAndAbsorbingPowers(CTE, CTM, 2*pi/lambda, b);
end

figure;
hold;
set(gca, 'fontsize', 16, 'fontname', 'times');
plot(lambda1/nm, PTEabs1./PTE0abs1, '-r', 'LineWidth', 2.5);
plot(lambda1/nm, PTMabs1./PTM0abs1, ':r', 'LineWidth', 2.5);
xlabel('\lambda_0 (nm)');
legend('TM', 'TE');
grid;


%-----------------------------------------------------------------ALUMINIUM
lambda0 = 495*nm;
b       = 0.167*lambda0;
a       = 0.95*b;

lambdamin = lambda0-200*nm;
lambdamax = lambda0+200*nm;
lambdapoints = 400;
lambda1   = lambdamin + (lambdamax-lambdamin)*[0:lambdapoints]/lambdapoints;

for lambdaindex = 0:lambdapoints
    lambda      = lambdamin + (lambdamax-lambdamin)*lambdaindex/lambdapoints;
    e1          = 1;
    e2          = AluminiumComplexPermittivity(c/lambda);
    [CTE, CTM]  = ScatteringCoefficients(e1, e2, a, b, 2*pi/lambda, N);
    [PTEabs1(1, lambdaindex+1),  PTEscat,...
     PTMabs1(1, lambdaindex+1),  PTMscat,...
    PTE0abs1(1, lambdaindex+1), PTE0scat,...
    PTM0abs1(1, lambdaindex+1), PTM0scat] =...
        NewScatteringAndAbsorbingPowers(CTE, CTM, 2*pi/lambda, b);
end

figure;
hold;
set(gca, 'fontsize', 16, 'fontname', 'times');
plot(lambda1/nm, PTEabs1./PTE0abs1, '-k', 'LineWidth', 2.5);
plot(lambda1/nm, PTMabs1./PTM0abs1, ':k', 'LineWidth', 2.5);
xlabel('\lambda_0 (nm)');
legend('TM', 'TE');
grid;