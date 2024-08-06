function [CTE, CTM] = ScatteringCoefficients(e1, e2, a, b, k0, N)
n = [-N:N];

k1 = k0.*sqrt(e1);
k2 = k0.*sqrt(e2);

I=1i;

CTE=(k1.*DJ(n,a.*k1).*(-(k0.*DJ(n,b.*k0).*H(n,b.*k2).*J(n,a.*k2))+k2.*J(n,b.*k0).*(-(DJ(n,b.*k2).*H(n,a.*k2))+DH(n,b.*k2).*J(n,a.*k2))+k0.*DJ(n,b.*k0).*H(n,a.*k2).*J(n,b.*k2))+k2.*J(n,a.*k1).*(k2.*(-(DH(n,b.*k2).*DJ(n,a.*k2))+DH(n,a.*k2).*DJ(n,b.*k2)).*J(n,b.*k0)+k0.*DJ(n,b.*k0).*(DJ(n,a.*k2).*H(n,b.*k2)-DH(n,a.*k2).*J(n,b.*k2))))./...
(I.^n.*(k2.*J(n,a.*k1).*(k2.*DH(n,b.*k2).*DJ(n,a.*k2).*H(n,b.*k0)-k0.*DH(n,b.*k0).*DJ(n,a.*k2).*H(n,b.*k2)+DH(n,a.*k2).*(-(k2.*DJ(n,b.*k2).*H(n,b.*k0))+k0.*DH(n,b.*k0).*J(n,b.*k2)))+k1.*DJ(n,a.*k1).*(k2.*DJ(n,b.*k2).*H(n,b.*k0).*H(n,a.*k2)-k2.*DH(n,b.*k2).*H(n,b.*k0).*J(n,a.*k2)+k0.*DH(n,b.*k0).*(H(n,b.*k2).*J(n,a.*k2)-H(n,a.*k2).*J(n,b.*k2)))));

CTM=(k2.*DJ(n,a.*k1).*(-(k2.*DJ(n,b.*k0).*H(n,b.*k2).*J(n,a.*k2))+k0.*J(n,b.*k0).*(-(DJ(n,b.*k2).*H(n,a.*k2))+DH(n,b.*k2).*J(n,a.*k2))+k2.*DJ(n,b.*k0).*H(n,a.*k2).*J(n,b.*k2))+k1.*J(n,a.*k1).*(k0.*(-(DH(n,b.*k2).*DJ(n,a.*k2))+DH(n,a.*k2).*DJ(n,b.*k2)).*J(n,b.*k0)+k2.*DJ(n,b.*k0).*(DJ(n,a.*k2).*H(n,b.*k2)-DH(n,a.*k2).*J(n,b.*k2))))./...
(I.^n.*(k1.*J(n,a.*k1).*(k0.*DH(n,b.*k2).*DJ(n,a.*k2).*H(n,b.*k0)-k2.*DH(n,b.*k0).*DJ(n,a.*k2).*H(n,b.*k2)+DH(n,a.*k2).*(-(k0.*DJ(n,b.*k2).*H(n,b.*k0))+k2.*DH(n,b.*k0).*J(n,b.*k2)))+k2.*DJ(n,a.*k1).*(k0.*DJ(n,b.*k2).*H(n,b.*k0).*H(n,a.*k2)-k0.*DH(n,b.*k2).*H(n,b.*k0).*J(n,a.*k2)+k2.*DH(n,b.*k0).*(H(n,b.*k2).*J(n,a.*k2)-H(n,a.*k2).*J(n,b.*k2)))));

function [out] = J(n, z)
out = besselj(n, z);

function [out] = H(n, z)
out = besselh(n, 2, z);

function [out] = DJ(n, z)
out = (1/2)*(J(-1+n,z)-J(1+n,z));

function [out] = DH(n, z)
out = (1/2)*(H(-1+n,z)-H(1+n,z));
