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
    Al_Air_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Al_Air_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Al_Air_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_Air_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Al_Air_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Al_Air_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_Air_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Al_Air_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Al_Air_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Al_Air_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Al_Air_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Al_Air_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Al_Air_Data.mat',...
    'Al_Air_PmaxTEabs',  'Al_Air_aoboptTEabs',  'Al_Air_boloptTEabs',...
    'Al_Air_PmaxTMabs',  'Al_Air_aoboptTMabs',  'Al_Air_boloptTMabs',...
    'Al_Air_PmaxTEscat', 'Al_Air_aoboptTEscat', 'Al_Air_boloptTEscat',...
    'Al_Air_PmaxTMscat', 'Al_Air_aoboptTMscat', 'Al_Air_boloptTMscat',...
    '-mat');

for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = AluminiumComplexPermittivity(f);
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
    Air_Al_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Air_Al_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Air_Al_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Air_Al_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Air_Al_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Air_Al_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Air_Al_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Air_Al_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Air_Al_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Air_Al_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Air_Al_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Air_Al_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Air_Al_Data.mat',...
    'Air_Al_PmaxTEabs',  'Air_Al_aoboptTEabs',  'Air_Al_boloptTEabs',...
    'Air_Al_PmaxTMabs',  'Air_Al_aoboptTMabs',  'Air_Al_boloptTMabs',...
    'Air_Al_PmaxTEscat', 'Air_Al_aoboptTEscat', 'Air_Al_boloptTEscat',...
    'Air_Al_PmaxTMscat', 'Air_Al_aoboptTMscat', 'Air_Al_boloptTMscat',...
    '-mat');