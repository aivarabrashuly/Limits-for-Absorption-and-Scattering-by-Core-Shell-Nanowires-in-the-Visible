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
    [e2, flag2] = GalliumAntimonideComplexPermittivity(f);
    
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
    Ag_GaSb_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Ag_GaSb_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ag_GaSb_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ag_GaSb_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Ag_GaSb_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ag_GaSb_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ag_GaSb_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Ag_GaSb_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ag_GaSb_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ag_GaSb_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Ag_GaSb_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ag_GaSb_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Ag_GaSb_Data.mat',...
    'Ag_GaSb_PmaxTEabs',  'Ag_GaSb_aoboptTEabs',  'Ag_GaSb_boloptTEabs',...
    'Ag_GaSb_PmaxTMabs',  'Ag_GaSb_aoboptTMabs',  'Ag_GaSb_boloptTMabs',...
    'Ag_GaSb_PmaxTEscat', 'Ag_GaSb_aoboptTEscat', 'Ag_GaSb_boloptTEscat',...
    'Ag_GaSb_PmaxTMscat', 'Ag_GaSb_aoboptTMscat', 'Ag_GaSb_boloptTMscat',...
    '-mat');

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = SilverComplexPermittivity(f);
    [e1, flag2] = GalliumAntimonideComplexPermittivity(f);
    
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
    GaSb_Ag_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    GaSb_Ag_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    GaSb_Ag_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaSb_Ag_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    GaSb_Ag_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    GaSb_Ag_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaSb_Ag_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    GaSb_Ag_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    GaSb_Ag_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaSb_Ag_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    GaSb_Ag_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    GaSb_Ag_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('GaSb_Ag_Data.mat',...
    'GaSb_Ag_PmaxTEabs',  'GaSb_Ag_aoboptTEabs',  'GaSb_Ag_boloptTEabs',...
    'GaSb_Ag_PmaxTMabs',  'GaSb_Ag_aoboptTMabs',  'GaSb_Ag_boloptTMabs',...
    'GaSb_Ag_PmaxTEscat', 'GaSb_Ag_aoboptTEscat', 'GaSb_Ag_boloptTEscat',...
    'GaSb_Ag_PmaxTMscat', 'GaSb_Ag_aoboptTMscat', 'GaSb_Ag_boloptTMscat',...
    '-mat');