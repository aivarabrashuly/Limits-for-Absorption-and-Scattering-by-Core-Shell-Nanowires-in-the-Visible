clear all;

nm      = 10^(-9);
c       = 3*10^8;
N       = 6;

lambdamin = 400*nm;
lambdamax = 700*nm;
lambdapoints = 300;
lambda1   = lambdamin + (lambdamax-lambdamin)*[0:lambdapoints]/lambdapoints;

%---------------------------------------TE (TM IN PAPER)---------AlSb ALONE
b       = 100*nm;
a       = 90*nm;

for lambdaindex = 0:lambdapoints
    lambda      = lambdamin + (lambdamax-lambdamin)*lambdaindex/lambdapoints;
    e1          = SilverComplexPermittivity(c/lambda);
    e2          = AmorphousSiliconComplexPermittivity(c/lambda);
    [CTE, CTM]  = ScatteringCoefficients(e1, e2, a, b, 2*pi/lambda, N);
    [PTEabs1(1, lambdaindex+1),  PTEscat1(1, lambdaindex+1),...
     PTMabs1(1, lambdaindex+1),  PTMscat1(1, lambdaindex+1),...
    PTE0abs1(1, lambdaindex+1), PTE0scat1(1, lambdaindex+1),...
    PTM0abs1(1, lambdaindex+1), PTM0scat1(1, lambdaindex+1)] =...
        NewScatteringAndAbsorbingPowers(CTE, CTM, 2*pi/lambda, b);
end

figure;
hold;
set(gca, 'fontsize', 16, 'fontname', 'times');
plot(lambda1/nm, (PTEabs1./PTE0abs1)*(2*b)/nm, '-r', 'LineWidth', 3.1);
plot(lambda1/nm, (PTMabs1./PTM0abs1)*(2*b)/nm, '-b', 'LineWidth', 3.1);
plot(lambda1/nm, (lambda1/(2*pi))/nm, '--k', 'LineWidth', 3.1);
xlabel('\lambda_0 (nm)');
legend('TM abs', 'TE abs');
grid;

figure;
hold;
set(gca, 'fontsize', 16, 'fontname', 'times');
plot(lambda1/nm, PTEabs1./PTE0abs1, '-r', 'LineWidth', 3.1);
plot(lambda1/nm, PTMabs1./PTM0abs1, '-b', 'LineWidth', 3.1);
% plot(lambda1/nm, (lambda1/(2*pi))/nm, '--k', 'LineWidth', 3.1);
xlabel('\lambda_0 (nm)');
legend('TM abs', 'TE abs');
grid;