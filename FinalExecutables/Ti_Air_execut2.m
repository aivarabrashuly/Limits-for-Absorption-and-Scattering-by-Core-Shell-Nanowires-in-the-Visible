clear all;

c       = 3*10^8;
nm      = 10^(-9);
N       = 12;

bolmin        = 0.001;
bolmax        = 0.25;
bolpoints     = 200;
bol1          = bolmin + (bolmax-bolmin)*[0:bolpoints]/bolpoints;

aobmin        = 0.001;
aobmax        = 0.999;
aobpoints     = 200;
aob1          = aobmin + (aobmax-aobmin)*[0:aobpoints]/aobpoints;

lambda0min    = 300*nm;
lambda0max    = 800*nm;
lambda0points = 50;
lambda01      = lambda0min + (lambda0max-lambda0min)*[0:lambda0points]/lambda0points;

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e1, flag1] = TitaniumComplexPermittivity(f);
    e2          = 1;
    
    for bolindex = 0:bolpoints
        bol = bolmin + (bolmax-bolmin)*bolindex/bolpoints;
        b   = bol*lambda0;
        for aobindex = 0:aobpoints
            aob = aobmin + (aobmax-aobmin)*aobindex/aobpoints;
            a   = aob*b;
            [CTE, CTM] = ScatteringCoefficients(e1, e2, a, b, k0, N);
            [PTEabs,  PTEscat,  PTMabs,  PTMscat,...
                PTE0abs, PTE0scat, PTM0abs, PTM0scat] =...
                NewScatteringAndAbsorbingPowers(CTE, CTM, k0, b);
            NormalizedPTEabs2(bolindex+1, aobindex+1) = PTEabs/PTE0abs;
            NormalizedPTMabs2(bolindex+1, aobindex+1) = PTMabs/PTM0abs;
            NormalizedPTEscat2(bolindex+1, aobindex+1) = PTEscat/PTE0scat;
            NormalizedPTMscat2(bolindex+1, aobindex+1) = PTMscat/PTM0scat;
        end
    end
    
    Q = NormalizedPTEabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_Air_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Ti_Air_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ti_Air_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_Air_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Ti_Air_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ti_Air_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_Air_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Ti_Air_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ti_Air_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_Air_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Ti_Air_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ti_Air_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Ti_Air_Data.mat',...
    'Ti_Air_PmaxTEabs',  'Ti_Air_aoboptTEabs',  'Ti_Air_boloptTEabs',...
    'Ti_Air_PmaxTMabs',  'Ti_Air_aoboptTMabs',  'Ti_Air_boloptTMabs',...
    'Ti_Air_PmaxTEscat', 'Ti_Air_aoboptTEscat', 'Ti_Air_boloptTEscat',...
    'Ti_Air_PmaxTMscat', 'Ti_Air_aoboptTMscat', 'Ti_Air_boloptTMscat',...
    '-mat');

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = TitaniumComplexPermittivity(f);
    e1          = 1;
    
    for bolindex = 0:bolpoints
        bol = bolmin + (bolmax-bolmin)*bolindex/bolpoints;
        b   = bol*lambda0;
        for aobindex = 0:aobpoints
            aob = aobmin + (aobmax-aobmin)*aobindex/aobpoints;
            a   = aob*b;
            [CTE, CTM] = ScatteringCoefficients(e1, e2, a, b, k0, N);
            [PTEabs,  PTEscat,  PTMabs,  PTMscat,...
                PTE0abs, PTE0scat, PTM0abs, PTM0scat] =...
                NewScatteringAndAbsorbingPowers(CTE, CTM, k0, b);
            NormalizedPTEabs2(bolindex+1, aobindex+1) = PTEabs/PTE0abs;
            NormalizedPTMabs2(bolindex+1, aobindex+1) = PTMabs/PTM0abs;
            NormalizedPTEscat2(bolindex+1, aobindex+1) = PTEscat/PTE0scat;
            NormalizedPTMscat2(bolindex+1, aobindex+1) = PTMscat/PTM0scat;
        end
    end
    
    Q = NormalizedPTEabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Air_Ti_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Air_Ti_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Air_Ti_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Air_Ti_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Air_Ti_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Air_Ti_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Air_Ti_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Air_Ti_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Air_Ti_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Air_Ti_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Air_Ti_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Air_Ti_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Air_Ti_Data.mat',...
    'Air_Ti_PmaxTEabs',  'Air_Ti_aoboptTEabs',  'Air_Ti_boloptTEabs',...
    'Air_Ti_PmaxTMabs',  'Air_Ti_aoboptTMabs',  'Air_Ti_boloptTMabs',...
    'Air_Ti_PmaxTEscat', 'Air_Ti_aoboptTEscat', 'Air_Ti_boloptTEscat',...
    'Air_Ti_PmaxTMscat', 'Air_Ti_aoboptTMscat', 'Air_Ti_boloptTMscat',...
    '-mat');