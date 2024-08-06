function [er, flag] = IndiumArsenideComplexPermittivity(f)

%D. E. Aspnes and A. A. Studna. Dielectric functions and optical parameters 
% of Si, Ge, GaP, GaAs, GaSb, InP, InAs, and InSb from 1.5 to 6.0 eV, 
% Phys. Rev. B 27, 985-1009 (1983)

flag = 0;

c = 3*10^8;

%The first column is wavelength in um!!!

Data = [0.2066 1.434 2.112
        0.2101 1.383 2.084
        0.2138 1.333 2.102
        0.2175 1.293 2.163
        0.2214 1.276 2.248
        0.2254 1.282 2.344
        0.2296 1.312 2.449
        0.2339 1.366 2.555
        0.2384 1.436 2.646
        0.2431 1.484 2.732
        0.2480 1.524 2.871
        0.2530 1.608 3.081
        0.2583 1.803 3.349
        0.2638 2.205 3.575
        0.2695 2.705 3.581
        0.2755 3.194 3.445
        0.2818 3.644 3.042
        0.2883 3.761 2.478
        0.2952 3.615 2.099
        0.3024 3.449 1.903
        0.3100 3.313 1.799
        0.3179 3.208 1.743
        0.3263 3.129 1.719
        0.3351 3.069 1.715
        0.3444 3.030 1.728
        0.3542 3.008 1.754
        0.3647 3.004 1.790
        0.3757 3.018 1.836
        0.3875 3.051 1.891
        0.3999 3.108 1.957
        0.4133 3.197 2.034
        0.4275 3.337 2.129
        0.4428 3.626 2.208
        0.4592 3.911 2.016
        0.4769 4.021 1.885
        0.4959 4.364 1.786
        0.5166 4.466 1.283
        0.5391 4.331 0.991
        0.5636 4.199 0.822
        0.5904 4.088 0.712
        0.6199 3.995 0.634
        0.6525 3.917 0.572
        0.6888 3.851 0.530
        0.7293 3.798 0.493
        0.7749 3.755 0.463
        0.8266 3.714 0.432];

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