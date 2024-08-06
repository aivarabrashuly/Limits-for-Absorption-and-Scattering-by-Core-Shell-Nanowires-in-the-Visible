clear all;

nm      = 10^(-9);
c       = 3*10^8;
N       = 6;
lambda0 = 700*nm;


lambdamin = lambda0-100*nm;
lambdamax = lambda0+100*nm;
lambdapoints = 200;
lambda1   = lambdamin + (lambdamax-lambdamin)*[0:lambdapoints]/lambdapoints;

%---------------------------------------TE (TM IN PAPER)---------AlSb ALONE
b       = 0.098*lambda0;
a       = 0.1*b;

for lambdaindex = 0:lambdapoints
    lambda      = lambdamin + (lambdamax-lambdamin)*lambdaindex/lambdapoints;
    e1          = AluminiumAntimonideComplexPermittivity(c/lambda);
    e2          = AluminiumAntimonideComplexPermittivity(c/lambda);
    [CTE, CTM]  = ScatteringCoefficients(e1, e2, a, b, 2*pi/lambda, N);
    [PTEabs1,  PTEscat1(1, lambdaindex+1),...
     PTMabs1,  PTMscat1(1, lambdaindex+1),...
    PTE0abs1, PTE0scat1(1, lambdaindex+1),...
    PTM0abs1, PTM0scat1(1, lambdaindex+1)] =...
        NewScatteringAndAbsorbingPowers(CTE, CTM, 2*pi/lambda, b);
end


figure;
hold;
set(gca, 'fontsize', 16, 'fontname', 'times');
% plot(lambda1/nm, PTEscat1./PTE0scat1, '-r', 'LineWidth', 2.1);
% plot(lambda1/nm, PTMscat1./PTM0scat1, '-b', 'LineWidth', 2.1);

%---------------------------------------------TE (TM IN PAPER)------Ag/AlSb
b       = 0.099*lambda0;
a       = 0.24*b;

for lambdaindex = 0:lambdapoints
    lambda      = lambdamin + (lambdamax-lambdamin)*lambdaindex/lambdapoints;
    e1          = SilverComplexPermittivity(c/lambda);
    e2          = AluminiumAntimonideComplexPermittivity(c/lambda);
    [CTE, CTM]  = ScatteringCoefficients(e1, e2, a, b, 2*pi/lambda, N);
    [PTEabs1,  PTEscat1(1, lambdaindex+1),...
     PTMabs1,  PTMscat1(1, lambdaindex+1),...
    PTE0abs1, PTE0scat1(1, lambdaindex+1),...
    PTM0abs1, PTM0scat1(1, lambdaindex+1)] =...
        NewScatteringAndAbsorbingPowers(CTE, CTM, 2*pi/lambda, b);
end

% figure;
% hold;
% set(gca, 'fontsize', 16, 'fontname', 'times');
% plot(lambda1/nm, PTEscat1./PTE0scat1, '-g', 'LineWidth', 2.1);
% plot(lambda1/nm, PTMscat1./PTM0scat1, '-k', 'LineWidth', 2.1);

%---------------------------------------------TM (TE IN PAPER)------Ag/AlSb
b       = 0.188*lambda0;
a       = 0.81*b;

for lambdaindex = 0:lambdapoints
    lambda      = lambdamin + (lambdamax-lambdamin)*lambdaindex/lambdapoints;
    e1          = SilverComplexPermittivity(c/lambda);
    e2          = AluminiumAntimonideComplexPermittivity(c/lambda);
    [CTE, CTM]  = ScatteringCoefficients(e1, e2, a, b, 2*pi/lambda, N);
    [PTEabs1,  PTEscat1(1, lambdaindex+1),...
     PTMabs1,  PTMscat1(1, lambdaindex+1),...
    PTE0abs1, PTE0scat1(1, lambdaindex+1),...
    PTM0abs1, PTM0scat1(1, lambdaindex+1)] =...
        NewScatteringAndAbsorbingPowers(CTE, CTM, 2*pi/lambda, b);
end

% figure;
% hold;
% set(gca, 'fontsize', 16, 'fontname', 'times');
% plot(lambda1/nm, PTEscat1./PTE0scat1, '--r', 'LineWidth', 2.1);
% plot(lambda1/nm, PTMscat1./PTM0scat1, '--b', 'LineWidth', 2.1);

%---------------------------------------------TM (TE IN PAPER)------Au/AlSb
b       = 0.194*lambda0;
a       = 0.84*b;

for lambdaindex = 0:lambdapoints
    lambda      = lambdamin + (lambdamax-lambdamin)*lambdaindex/lambdapoints;
    e1          = GoldComplexPermittivity(c/lambda);
    e2          = AluminiumAntimonideComplexPermittivity(c/lambda);
    [CTE, CTM]  = ScatteringCoefficients(e1, e2, a, b, 2*pi/lambda, N);
    [PTEabs1,  PTEscat1(1, lambdaindex+1),...
     PTMabs1,  PTMscat1(1, lambdaindex+1),...
    PTE0abs1, PTE0scat1(1, lambdaindex+1),...
    PTM0abs1, PTM0scat1(1, lambdaindex+1)] =...
        NewScatteringAndAbsorbingPowers(CTE, CTM, 2*pi/lambda, b);
end

% figure;
% hold;
% set(gca, 'fontsize', 16, 'fontname', 'times');
plot(lambda1/nm, PTEscat1./PTE0scat1, '--g', 'LineWidth', 2.1);
plot(lambda1/nm, PTMscat1./PTM0scat1, '--k', 'LineWidth', 2.1);