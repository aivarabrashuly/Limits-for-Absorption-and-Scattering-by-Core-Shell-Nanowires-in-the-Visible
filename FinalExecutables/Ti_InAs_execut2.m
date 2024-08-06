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
    [e2, flag2] = IndiumArsenideComplexPermittivity(f);
    
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
    Ti_InAs_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Ti_InAs_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ti_InAs_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_InAs_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Ti_InAs_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ti_InAs_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_InAs_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Ti_InAs_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ti_InAs_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_InAs_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Ti_InAs_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ti_InAs_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Ti_InAs_Data.mat',...
    'Ti_InAs_PmaxTEabs',  'Ti_InAs_aoboptTEabs',  'Ti_InAs_boloptTEabs',...
    'Ti_InAs_PmaxTMabs',  'Ti_InAs_aoboptTMabs',  'Ti_InAs_boloptTMabs',...
    'Ti_InAs_PmaxTEscat', 'Ti_InAs_aoboptTEscat', 'Ti_InAs_boloptTEscat',...
    'Ti_InAs_PmaxTMscat', 'Ti_InAs_aoboptTMscat', 'Ti_InAs_boloptTMscat',...
    '-mat');

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = TitaniumComplexPermittivity(f);
    [e1, flag2] = IndiumArsenideComplexPermittivity(f);
    
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
    InAs_Ti_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    InAs_Ti_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    InAs_Ti_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    InAs_Ti_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    InAs_Ti_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    InAs_Ti_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    InAs_Ti_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    InAs_Ti_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    InAs_Ti_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    InAs_Ti_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    InAs_Ti_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    InAs_Ti_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('InAs_Ti_Data.mat',...
    'InAs_Ti_PmaxTEabs',  'InAs_Ti_aoboptTEabs',  'InAs_Ti_boloptTEabs',...
    'InAs_Ti_PmaxTMabs',  'InAs_Ti_aoboptTMabs',  'InAs_Ti_boloptTMabs',...
    'InAs_Ti_PmaxTEscat', 'InAs_Ti_aoboptTEscat', 'InAs_Ti_boloptTEscat',...
    'InAs_Ti_PmaxTMscat', 'InAs_Ti_aoboptTMscat', 'InAs_Ti_boloptTMscat',...
    '-mat');