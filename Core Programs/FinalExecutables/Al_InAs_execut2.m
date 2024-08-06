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
    Al_InAs_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Al_InAs_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Al_InAs_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_InAs_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Al_InAs_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Al_InAs_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_InAs_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Al_InAs_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Al_InAs_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_InAs_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Al_InAs_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Al_InAs_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Al_InAs_Data.mat',...
    'Al_InAs_PmaxTEabs',  'Al_InAs_aoboptTEabs',  'Al_InAs_boloptTEabs',...
    'Al_InAs_PmaxTMabs',  'Al_InAs_aoboptTMabs',  'Al_InAs_boloptTMabs',...
    'Al_InAs_PmaxTEscat', 'Al_InAs_aoboptTEscat', 'Al_InAs_boloptTEscat',...
    'Al_InAs_PmaxTMscat', 'Al_InAs_aoboptTMscat', 'Al_InAs_boloptTMscat',...
    '-mat');

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = AluminiumComplexPermittivity(f);
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
    InAs_Al_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    InAs_Al_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    InAs_Al_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    InAs_Al_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    InAs_Al_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    InAs_Al_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    InAs_Al_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    InAs_Al_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    InAs_Al_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    InAs_Al_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    InAs_Al_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    InAs_Al_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('InAs_Al_Data.mat',...
    'InAs_Al_PmaxTEabs',  'InAs_Al_aoboptTEabs',  'InAs_Al_boloptTEabs',...
    'InAs_Al_PmaxTMabs',  'InAs_Al_aoboptTMabs',  'InAs_Al_boloptTMabs',...
    'InAs_Al_PmaxTEscat', 'InAs_Al_aoboptTEscat', 'InAs_Al_boloptTEscat',...
    'InAs_Al_PmaxTMscat', 'InAs_Al_aoboptTMscat', 'InAs_Al_boloptTMscat',...
    '-mat');