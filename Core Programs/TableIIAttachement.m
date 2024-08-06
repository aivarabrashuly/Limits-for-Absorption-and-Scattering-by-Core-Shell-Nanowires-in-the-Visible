clear all;

nm      = 10^(-9);
c       = 3*10^8;

%---------------------------------------------------------------DIPOLE-LIKE 
lambda0 = 620*nm;
k0      = 2*pi/lambda0;
f       = c/lambda0;
b       = 0.081*lambda0;
a       = 0.49*b;
N       = 6;
e1      = SilverComplexPermittivity(f);
e2      = AluminiumAntimonideComplexPermittivity(f);

[CTE, CTM] = ScatteringCoefficients(e1, e2, a, b, k0, N);

S0 = CTM(1, N+1);
S1 = CTM(1, N+2);
ReS0 = real(S0);
ImS0 = imag(S0);
ReS1 = real(S1);
ImS1 = imag(S1);

disp(abs(S0)^2);
disp(4*abs(S1)^2);
disp(4*(ImS0*ReS1-ImS1*ReS0));
disp(angle(S0));
disp(angle(S1));

figure;
hold;
plot([-N:N], real(CTM), '-b', 'LineWidth', 2.8);
plot([-N:N], imag(CTM), '-r', 'LineWidth', 2.8);
plot([-N:N], abs(CTM), '-k',  'LineWidth', 2.8);
xlabel('n');
ylabel('|S_n^{TE}|');
title('Dipole-Like Pattern');
grid;

[PTEabs,  PTEscat,  PTMabs,  PTMscat,...
          PTE0abs, PTE0scat, PTM0abs, PTM0scat] =...
          NewScatteringAndAbsorbingPowers(CTE, CTM, k0, b);
      
 disp(PTMabs);
 disp(PTM0abs);

% %--------------------------------------------------------------HUYGENS-LIKE 
% lambda0 = 500*nm;
% k0      = 2*pi/lambda0;
% f       = c/lambda0;
% b       = 0.138*lambda0;
% a       = 0.68*b;
% N       = 6;
% e1      = TitaniumComplexPermittivity(f);
% e2      = GalliumPhosphideComplexPermittivity(f);
% 
% [CTE, CTM] = ScatteringCoefficients(e1, e2, a, b, k0, N);
% 
% S0 = CTM(1, N+1);
% S1 = CTM(1, N+2);
% ReS0 = real(S0);
% ImS0 = imag(S0);
% ReS1 = real(S1);
% ImS1 = imag(S1);
% 
% % disp(abs(S0)^2);
% % disp(4*abs(S1)^2);
% % disp(4*(ImS0*ReS1-ImS1*ReS0));
% % disp(angle(S0));
% % disp(angle(S1));
% 
% 
% [PTEabs,  PTEscat,  PTMabs,  PTMscat,...
%           PTE0abs, PTE0scat, PTM0abs, PTM0scat] =...
%           NewScatteringAndAbsorbingPowers(CTE, CTM, k0, b);
% 
% figure;
% hold;
% plot([-N:N], real(CTM), '-b', 'LineWidth', 2.8);
% plot([-N:N], imag(CTM), '-r', 'LineWidth', 2.8);
% plot([-N:N], abs(CTM), '-k',  'LineWidth', 2.8);
% xlabel('n');
% ylabel('|S_n^{TE}|');
% title('Huygens-Like Pattern');
% grid;


clear all;

nm      = 10^(-9);
c       = 3*10^8;

% %---------------------------------------------------------------SILVER TUBE 
% lambda0 = 620*nm;
% k0      = 2*pi/lambda0;
% f       = c/lambda0;
% b       = 0.081*lambda0;
% a       = 0.49*b;
% N       = 6;
% e1      = SilverComplexPermittivity(f);
% e2      = AluminiumAntimonideComplexPermittivity(f);
% 
% [CTE, CTM] = ScatteringCoefficients(e1, e2, a, b, k0, N);
% 
% S0 = CTM(1, N+1);
% S1 = CTM(1, N+2);
% ReS0 = real(S0);
% ImS0 = imag(S0);
% ReS1 = real(S1);
% ImS1 = imag(S1);
% 
% % disp(abs(S0)^2);
% % disp(4*abs(S1)^2);
% % disp(4*(ImS0*ReS1-ImS1*ReS0));
% % disp(angle(S0));
% % disp(angle(S1));
% 
% figure;
% hold;
% plot([-N:N], real(CTM), '-b', 'LineWidth', 2.8);
% plot([-N:N], imag(CTM), '-r', 'LineWidth', 2.8);
% plot([-N:N], abs(CTM), '-k',  'LineWidth', 2.8);
% xlabel('n');
% ylabel('|S_n^{TE}|');
% title('Dipole-Like Pattern');
% grid;
% 
% [PTEabs,  PTEscat,  PTMabs,  PTMscat,...
%           PTE0abs, PTE0scat, PTM0abs, PTM0scat] =...
%           NewScatteringAndAbsorbingPowers(CTE, CTM, k0, b);
%       
%  disp(PTMabs);
%  disp(PTM0abs);