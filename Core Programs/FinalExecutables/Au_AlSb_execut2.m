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
    [e1, flag1] = GoldComplexPermittivity(f);
    [e2, flag2] = AluminiumAntimonideComplexPermittivity(f);
    
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
    Au_AlSb_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    Au_AlSb_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Au_AlSb_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Au_AlSb_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    Au_AlSb_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    Au_AlSb_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Au_AlSb_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    Au_AlSb_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Au_AlSb_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    Au_AlSb_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    Au_AlSb_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    Au_AlSb_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('Au_AlSb_Data.mat',...
    'Au_AlSb_PmaxTEabs',  'Au_AlSb_aoboptTEabs',  'Au_AlSb_boloptTEabs',...
    'Au_AlSb_PmaxTMabs',  'Au_AlSb_aoboptTMabs',  'Au_AlSb_boloptTMabs',...
    'Au_AlSb_PmaxTEscat', 'Au_AlSb_aoboptTEscat', 'Au_AlSb_boloptTEscat',...
    'Au_AlSb_PmaxTMscat', 'Au_AlSb_aoboptTMscat', 'Au_AlSb_boloptTMscat',...
    '-mat');


for lambda0index = 0:lambda0points
    lambda0     = lambda0min + (lambda0max-lambda0min)*lambda0index/lambda0points;
    f           = c/lambda0;
    k0          = 2*pi/lambda0;
    [e2, flag1] = GoldComplexPermittivity(f);
    [e1, flag2] = AluminiumAntimonideComplexPermittivity(f);
    
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
    AlSb_Au_PmaxTEabs(lambda0index+1, 1)    = Qmax;
    AlSb_Au_aoboptTEabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    AlSb_Au_boloptTEabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTMabs2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    AlSb_Au_PmaxTMabs(lambda0index+1, 1)    = Qmax;
    AlSb_Au_aoboptTMabs(lambda0index+1, 1)  = aob1(1, aoboptindex);
    AlSb_Au_boloptTMabs(lambda0index+1, 1)  = bol1(1, boloptindex);
    
    Q = NormalizedPTEscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    AlSb_Au_PmaxTEscat(lambda0index+1, 1)   = Qmax;
    AlSb_Au_aoboptTEscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    AlSb_Au_boloptTEscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    Q = NormalizedPTMscat2;
    Qmax = MyMax(Q);
    [boloptindex, aoboptindex] = find(Q==Qmax);
    AlSb_Au_PmaxTMscat(lambda0index+1, 1)   = Qmax;
    AlSb_Au_aoboptTMscat(lambda0index+1, 1) = aob1(1, aoboptindex);
    AlSb_Au_boloptTMscat(lambda0index+1, 1) = bol1(1, boloptindex);
    
    disp(lambda0index);
end

save('AlSb_Au_Data.mat',...
    'AlSb_Au_PmaxTEabs',  'AlSb_Au_aoboptTEabs',  'AlSb_Au_boloptTEabs',...
    'AlSb_Au_PmaxTMabs',  'AlSb_Au_aoboptTMabs',  'AlSb_Au_boloptTMabs',...
    'AlSb_Au_PmaxTEscat', 'AlSb_Au_aoboptTEscat', 'AlSb_Au_boloptTEscat',...
    'AlSb_Au_PmaxTMscat', 'AlSb_Au_aoboptTMscat', 'AlSb_Au_boloptTMscat',...
    '-mat');