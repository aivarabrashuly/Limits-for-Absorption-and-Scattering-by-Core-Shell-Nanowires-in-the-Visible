function [er, flag] = TitaniumComplexPermittivity(f)

% P. B. Johnson and R. W. Christy. Optical constants of transition metals:
% Ti, V, Cr, Mn, Fe, Co, Ni, and Pd, Phys. Rev. B 9, 5056-5070 (1974)

flag = 0;

c = 3*10^8;
%The first column is wavelength in um!!!

Data = [0.188 1.10 1.62
        0.192 1.16 1.64
        0.195 1.22 1.66
        0.199 1.25 1.68
        0.203 1.27 1.69
        0.207 1.31 1.69
        0.212 1.31 1.68
        0.216 1.32 1.67
        0.221 1.32 1.66
        0.226 1.32 1.66
        0.231 1.31 1.68
        0.237 1.30 1.72
        0.243 1.28 1.77
        0.249 1.27 1.83
        0.255 1.26 1.91
        0.262 1.27 1.99
        0.269 1.27 2.07
        0.276 1.30 2.17
        0.284 1.35 2.26
        0.292 1.40 2.36
        0.301 1.45 2.46
        0.311 1.50 2.57
        0.320 1.55 2.66
        0.332 1.61 2.74
        0.342 1.72 2.82
        0.354 1.82 2.87
        0.368 1.90 2.90
        0.381 1.99 2.93
        0.397 2.08 2.95
        0.413 2.14 2.98
        0.431 2.21 3.01
        0.451 2.27 3.04
        0.471 2.32 3.10
        0.496 2.36 3.19
        0.521 2.44 3.30
        0.549 2.54 3.43
        0.582 2.60 3.58
        0.617 2.67 3.72
        0.659 2.76 3.84
        0.704 2.86 3.96
        0.756 3.00 4.01
        0.821 3.21 4.01
        0.892 3.29 3.96
        0.984 3.35 3.97
        1.088 3.50 4.02
        1.216 3.62 4.15
        1.393 3.67 4.37
        1.610 3.69 4.70
        1.937 3.51 5.19];

[N, UseLess] =  size(Data);
fmin   = c/(Data(N, 1)*10^(-6));
fmax   = c/(Data(1, 1)*10^(-6));

if (f>fmin)&&(f<fmax)
    v             = 3*10^8;
    lambda0       = v/f;
    [UseLess, l0indexa] = min(abs(Data(:, 1)-lambda0*10^6));
    la            = Data(l0indexa, 1);
    na            = Data(l0indexa, 2);
    ka            = Data(l0indexa, 3);
    era           = (na-1i*ka)^2;
    Data(l0indexa, 1) = 10^10;
    [UseLess, l0indexb] = min(abs(Data(:, 1)-lambda0*10^6));
    lb            = Data(l0indexb, 1);
    nb            = Data(l0indexb, 2);
    kb            = Data(l0indexb, 3);
    erb           = (nb-1i*kb)^2;
    er            = era + (erb-era)*(lambda0*10^6-la)/(lb-la);
    flag          = 1;
else
    flag          = 0;
    er            = 0;
end