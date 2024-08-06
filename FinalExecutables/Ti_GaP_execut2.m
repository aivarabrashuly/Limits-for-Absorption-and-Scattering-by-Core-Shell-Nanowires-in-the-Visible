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
    [e2, flag2] = GalliumPhosphideComplexPermittivity(f);
    
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
    Ti_GaP_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Ti_GaP_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ti_GaP_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_GaP_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Ti_GaP_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ti_GaP_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_GaP_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Ti_GaP_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ti_GaP_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ti_GaP_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Ti_GaP_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ti_GaP_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Ti_GaP_Data.mat',...
    'Ti_GaP_PmaxTEabs',  'Ti_GaP_aoboptTEabs',  'Ti_GaP_boloptTEabs',...
    'Ti_GaP_PmaxTMabs',  'Ti_GaP_aoboptTMabs',  'Ti_GaP_boloptTMabs',...
    'Ti_GaP_PmaxTEscat', 'Ti_GaP_aoboptTEscat', 'Ti_GaP_boloptTEscat',...
    'Ti_GaP_PmaxTMscat', 'Ti_GaP_aoboptTMscat', 'Ti_GaP_boloptTMscat',...
    '-mat');

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = TitaniumComplexPermittivity(f);
    [e1, flag2] = GalliumPhosphideComplexPermittivity(f);
    
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
    GaP_Ti_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    GaP_Ti_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    GaP_Ti_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaP_Ti_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    GaP_Ti_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    GaP_Ti_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaP_Ti_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    GaP_Ti_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    GaP_Ti_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaP_Ti_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    GaP_Ti_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    GaP_Ti_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('GaP_Ti_Data.mat',...
    'GaP_Ti_PmaxTEabs',  'GaP_Ti_aoboptTEabs',  'GaP_Ti_boloptTEabs',...
    'GaP_Ti_PmaxTMabs',  'GaP_Ti_aoboptTMabs',  'GaP_Ti_boloptTMabs',...
    'GaP_Ti_PmaxTEscat', 'GaP_Ti_aoboptTEscat', 'GaP_Ti_boloptTEscat',...
    'GaP_Ti_PmaxTMscat', 'GaP_Ti_aoboptTMscat', 'GaP_Ti_boloptTMscat',...
    '-mat');