function [er, flag] = AmorphousSiliconComplexPermittivity(f)

%1) D. T. Pierce and W. E. Spicer, Electronic structure of amorphous Si from 
%photoemission and optical studies, Phys. Rev. B 5, 3017-3029 (1972)
%2) Handbook of Optical Constants of Solids, Edward D. Palik, ed. Academic Press, Boston, 1985
%(ref. 2 provides numerical values for the graphical data reported in ref. 1)

flag = 0;

c = 3*10^8;

%The first column is wavelength in um!!!

Data = [ 1.033E-1 3.27E-1 7.26E-1
        1.078E-1 3.63E-1 8.47E-1
        1.127E-1 3.92E-1 9.46E-1
        1.181E-1 4.23E-1 1.04E+0
        1.240E-1 4.59E-1 1.14E+0
        1.305E-1 4.97E-1 1.24E+0
        1.378E-1 5.43E-1 1.35E+0
        1.459E-1 5.97E-1 1.47E+0
        1.550E-1 6.60E-1 1.60E+0
        1.653E-1 7.35E-1 1.74E+0
        1.771E-1 8.32E-1 1.89E+0
        1.907E-1 9.51E-1 2.07E+0
        2.066E-1 1.11E+0 2.28E+0
        2.254E-1 1.35E+0 2.51E+0
        2.480E-1 1.69E+0 2.76E+0
        2.583E-1 1.86E+0 2.85E+0
        2.695E-1 2.07E+0 2.93E+0
        2.818E-1 2.30E+0 2.99E+0
        2.952E-1 2.56E+0 3.04E+0
        3.100E-1 2.87E+0 3.06E+0
        3.263E-1 3.21E+0 3.00E+0
        3.444E-1 3.55E+0 2.88E+0
        3.543E-1 3.73E+0 2.79E+0
        3.647E-1 3.90E+0 2.66E+0
        3.875E-1 4.17E+0 2.38E+0
        4.133E-1 4.38E+0 2.02E+0
        4.428E-1 4.47E+0 1.64E+0
        4.769E-1 4.49E+0 1.28E+0
        4.960E-1 4.47E+0 1.12E+0
        5.166E-1 4.46E+0 9.69E-1
        5.636E-1 4.36E+0 6.90E-1
        6.199E-1 4.23E+0 4.61E-1
        6.526E-1 4.17E+0 3.63E-1
        6.888E-1 4.09E+0 2.71E-1
        7.293E-1 4.01E+0 1.99E-1
        7.749E-1 3.93E+0 1.36E-1
        8.266E-1 3.86E+0 8.12E-2
        8.856E-1 3.77E+0 4.01E-2
        9.538E-1 3.68E+0 0.00E+0
        1.033E+0 3.61E+0 0.00E+0
        1.127E+0 3.57E+0 0.00E+0
        1.240E+0 3.54E+0 0.00E+0
        1.378E+0 3.50E+0 0.00E+0
        1.550E+0 3.48E+0 0.00E+0
        1.771E+0 3.45E+0 0.00E+0
        2.066E+0 3.44E+0 0.00E+0];

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