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
    [e1, flag1] = SilverComplexPermittivity(f);
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
    Ag_GaP_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Ag_GaP_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ag_GaP_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ag_GaP_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Ag_GaP_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ag_GaP_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ag_GaP_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Ag_GaP_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ag_GaP_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ag_GaP_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Ag_GaP_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ag_GaP_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Ag_GaP_Data.mat',...
    'Ag_GaP_PmaxTEabs',  'Ag_GaP_aoboptTEabs',  'Ag_GaP_boloptTEabs',...
    'Ag_GaP_PmaxTMabs',  'Ag_GaP_aoboptTMabs',  'Ag_GaP_boloptTMabs',...
    'Ag_GaP_PmaxTEscat', 'Ag_GaP_aoboptTEscat', 'Ag_GaP_boloptTEscat',...
    'Ag_GaP_PmaxTMscat', 'Ag_GaP_aoboptTMscat', 'Ag_GaP_boloptTMscat',...
    '-mat');

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = SilverComplexPermittivity(f);
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
    GaP_Ag_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    GaP_Ag_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    GaP_Ag_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaP_Ag_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    GaP_Ag_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    GaP_Ag_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaP_Ag_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    GaP_Ag_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    GaP_Ag_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaP_Ag_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    GaP_Ag_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    GaP_Ag_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('GaP_Ag_Data.mat',...
    'GaP_Ag_PmaxTEabs',  'GaP_Ag_aoboptTEabs',  'GaP_Ag_boloptTEabs',...
    'GaP_Ag_PmaxTMabs',  'GaP_Ag_aoboptTMabs',  'GaP_Ag_boloptTMabs',...
    'GaP_Ag_PmaxTEscat', 'GaP_Ag_aoboptTEscat', 'GaP_Ag_boloptTEscat',...
    'GaP_Ag_PmaxTMscat', 'GaP_Ag_aoboptTMscat', 'GaP_Ag_boloptTMscat',...
    '-mat');