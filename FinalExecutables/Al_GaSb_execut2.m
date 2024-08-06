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
    [e1, flag1] = AluminiumComplexPermittivity(f);
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
    Al_GaSb_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Al_GaSb_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Al_GaSb_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_GaSb_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Al_GaSb_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Al_GaSb_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_GaSb_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Al_GaSb_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Al_GaSb_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_GaSb_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Al_GaSb_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Al_GaSb_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Al_GaSb_Data.mat',...
    'Al_GaSb_PmaxTEabs',  'Al_GaSb_aoboptTEabs',  'Al_GaSb_boloptTEabs',...
    'Al_GaSb_PmaxTMabs',  'Al_GaSb_aoboptTMabs',  'Al_GaSb_boloptTMabs',...
    'Al_GaSb_PmaxTEscat', 'Al_GaSb_aoboptTEscat', 'Al_GaSb_boloptTEscat',...
    'Al_GaSb_PmaxTMscat', 'Al_GaSb_aoboptTMscat', 'Al_GaSb_boloptTMscat',...
    '-mat');

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = AluminiumComplexPermittivity(f);
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
    GaSb_Al_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    GaSb_Al_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    GaSb_Al_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaSb_Al_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    GaSb_Al_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    GaSb_Al_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaSb_Al_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    GaSb_Al_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    GaSb_Al_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    GaSb_Al_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    GaSb_Al_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    GaSb_Al_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('GaSb_Al_Data.mat',...
    'GaSb_Al_PmaxTEabs',  'GaSb_Al_aoboptTEabs',  'GaSb_Al_boloptTEabs',...
    'GaSb_Al_PmaxTMabs',  'GaSb_Al_aoboptTMabs',  'GaSb_Al_boloptTMabs',...
    'GaSb_Al_PmaxTEscat', 'GaSb_Al_aoboptTEscat', 'GaSb_Al_boloptTEscat',...
    'GaSb_Al_PmaxTMscat', 'GaSb_Al_aoboptTMscat', 'GaSb_Al_boloptTMscat',...
    '-mat');