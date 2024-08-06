function [PTEabs,  PTEscat,  PTMabs,  PTMscat,...
          PTE0abs, PTE0scat, PTM0abs, PTM0scat] =...
          NewScatteringAndAbsorbingPowers(CTE, CTM, k0, b)
n0=120*pi;
N=(length(CTE)-1)/2;
n=[-N:N];
ReCTE=real(CTE);
ReCTM=real(CTM);
ImCTE=imag(CTE);
ImCTM=imag(CTM);

CTE0 = besselj(n, k0*b)./besselh(n, 2, k0*b);
CTM0 = besselj(n, k0*b)./besselh(n, 2, k0*b);

PTEabs   = (-4/(k0*n0))*(sum(abs(CTE).^2)+sum(ReCTE.*cos(n*pi/2)-ImCTE.*sin(n*pi/2)));
PTEscat  = (+4/(k0*n0))*(sum(abs(CTE).^2));

PTMabs   = (-4*n0/k0)*(sum(abs(CTM).^2)+sum(ReCTM.*cos(n*pi/2)-ImCTM.*sin(n*pi/2)));
PTMscat  = (+4*n0/k0)*(sum(abs(CTM).^2));

PTE0abs  = 2*b/n0;
PTE0scat = (+4/(k0*n0))*(sum(abs(CTE0).^2));
PTM0abs  = 2*b*n0;
PTM0scat = (+4*n0/k0)*(sum(abs(CTM0).^2));