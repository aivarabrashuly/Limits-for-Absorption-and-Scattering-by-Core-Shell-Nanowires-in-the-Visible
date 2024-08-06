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
    [e2, flag2] = GermaniumComplexPermittivity(f);
    
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
    Al_Ge_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Al_Ge_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Al_Ge_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_Ge_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Al_Ge_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Al_Ge_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_Ge_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Al_Ge_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Al_Ge_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_Ge_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Al_Ge_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Al_Ge_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Al_Ge_Data.mat',...
    'Al_Ge_PmaxTEabs',  'Al_Ge_aoboptTEabs',  'Al_Ge_boloptTEabs',...
    'Al_Ge_PmaxTMabs',  'Al_Ge_aoboptTMabs',  'Al_Ge_boloptTMabs',...
    'Al_Ge_PmaxTEscat', 'Al_Ge_aoboptTEscat', 'Al_Ge_boloptTEscat',...
    'Al_Ge_PmaxTMscat', 'Al_Ge_aoboptTMscat', 'Al_Ge_boloptTMscat',...
    '-mat');

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = AluminiumComplexPermittivity(f);
    [e1, flag2] = GermaniumComplexPermittivity(f);
    
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
    Ge_Al_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Ge_Al_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ge_Al_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ge_Al_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Ge_Al_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Ge_Al_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ge_Al_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Ge_Al_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ge_Al_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Ge_Al_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Ge_Al_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Ge_Al_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Ge_Al_Data.mat',...
    'Ge_Al_PmaxTEabs',  'Ge_Al_aoboptTEabs',  'Ge_Al_boloptTEabs',...
    'Ge_Al_PmaxTMabs',  'Ge_Al_aoboptTMabs',  'Ge_Al_boloptTMabs',...
    'Ge_Al_PmaxTEscat', 'Ge_Al_aoboptTEscat', 'Ge_Al_boloptTEscat',...
    'Ge_Al_PmaxTMscat', 'Ge_Al_aoboptTMscat', 'Ge_Al_boloptTMscat',...
    '-mat');