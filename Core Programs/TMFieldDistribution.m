function [Hz2, x1, y1] =...
    TMFieldDistribution(xmin, xmax, ymin, ymax, a, b, e1, e2, lambda0)

N = 10;

k0 = 2*pi/lambda0;
k1 = k0*sqrt(e1);
k2 = k0*sqrt(e2);

PPWL = 600;

xpoints = round(PPWL*(xmax-xmin)/lambda0);
x1      = xmin + (xmax-xmin)*[0:xpoints]/xpoints;

ypoints = round(PPWL*(ymax-ymin)/lambda0);
y1      = ymin + (ymax-ymin)*[0:ypoints]/ypoints;

n1 = [-N:N];
Q1 = 1i.^(-n1);

[C11, C21, D21, C01] = ComputeCoefficients(n1, a, b, e1, e2, lambda0);

% figure; plot(abs(C01))

x2 = repmat(reshape(x1, [1, xpoints+1]), [ypoints+1, 1]);
y2 = repmat(reshape(y1, [ypoints+1, 1]), [1, xpoints+1]);

r2   = sqrt(x2.^2+y2.^2);
phi2 = atan2(y2, x2);

r3   = repmat(r2,   [1, 1, 2*N+1]);
phi3 = repmat(phi2, [1, 1, 2*N+1]);
n3   = repmat(reshape(n1, [1, 1, 2*N+1]),  [ypoints+1, xpoints+1, 1]);
C13  = repmat(reshape(C11, [1, 1, 2*N+1]), [ypoints+1, xpoints+1, 1]);
C23  = repmat(reshape(C21, [1, 1, 2*N+1]), [ypoints+1, xpoints+1, 1]);
D23  = repmat(reshape(D21, [1, 1, 2*N+1]), [ypoints+1, xpoints+1, 1]);
C03  = repmat(reshape(C01, [1, 1, 2*N+1]), [ypoints+1, xpoints+1, 1]);
Q3   = repmat(reshape(Q1, [1, 1, 2*N+1]),  [ypoints+1, xpoints+1, 1]);

H1z3 =  C13.*J(n3, k1.*r3).*exp(1i.*n3.*phi3);
H2z3 = (C23.*J(n3, k2.*r3)+D23.*H(n3, k2.*r3+eps)).*exp(1i.*n3.*phi3);
H0z3 = (Q3.*J(n3, k0.*r3)+C03.*H(n3, k0.*r3+eps)).*exp(1i.*n3.*phi3);

H1z2 = sum(H1z3, 3);
H2z2 = sum(H2z3, 3);
H0z2 = sum(H0z3, 3);

Hz2 =  H1z2.*(r2<=a) + H2z2.*(r2>a).*(r2<b) + H0z2.*(r2>=b);

function [C11, C21, D21, C01] = ComputeCoefficients(n, a, b, e1, e2, lambda0)

k0 = 2*pi/lambda0;
k1 = k0*sqrt(e1);
k2 = k0*sqrt(e2);

Q1 = 1i.^(-n);

I = 1i;

C11=(-((k2.*(DJ(n,b.*k0).*H(n,b.*k0)-DH(n,b.*k0).*J(n,b.*k0)).*(k2.*DJ(n,a.*k1).*H(n,a.*k2)-k1.*DH(n,a.*k2).*J(n,a.*k1)).*J(n,a.*k2))./(I.^n.*(-(k0.*k2.*DJ(n,a.*k1).*DJ(n,b.*k2).*H(n,b.*k0).*H(n,a.*k2))-k0.*k1.*DH(n,b.*k2).*DJ(n,a.*k2).*H(n,b.*k0).*J(n,a.*k1)+k0.*k1.*DH(n,a.*k2).*DJ(n,b.*k2).*H(n,b.*k0).*J(n,a.*k1)+k1.*k2.*DH(n,b.*k0).*DJ(n,a.*k2).*H(n,b.*k2).*J(n,a.*k1)+k0.*k2.*DH(n,b.*k2).*DJ(n,a.*k1).*H(n,b.*k0).*J(n,a.*k2)-k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,b.*k2).*J(n,a.*k2)+k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,a.*k2).*J(n,b.*k2)-k1.*k2.*DH(n,b.*k0).*DH(n,a.*k2).*J(n,a.*k1).*J(n,b.*k2))))+...
    (k2.*H(n,a.*k2).*(DJ(n,b.*k0).*H(n,b.*k0)-DH(n,b.*k0).*J(n,b.*k0)).*(k2.*DJ(n,a.*k1).*H(n,a.*k2)-k1.*DH(n,a.*k2).*J(n,a.*k1)).*(k1.*DJ(n,a.*k2).*J(n,a.*k1)-k2.*DJ(n,a.*k1).*J(n,a.*k2)))./(I.^n.*(-(k2.*DJ(n,a.*k1).*H(n,a.*k2))+k1.*DH(n,a.*k2).*J(n,a.*k1)).*(-(k0.*k2.*DJ(n,a.*k1).*DJ(n,b.*k2).*H(n,b.*k0).*H(n,a.*k2))-k0.*k1.*DH(n,b.*k2).*DJ(n,a.*k2).*H(n,b.*k0).*J(n,a.*k1)+k0.*k1.*DH(n,a.*k2).*DJ(n,b.*k2).*H(n,b.*k0).*J(n,a.*k1)+k1.*k2.*DH(n,b.*k0).*DJ(n,a.*k2).*H(n,b.*k2).*J(n,a.*k1)+k0.*k2.*DH(n,b.*k2).*DJ(n,a.*k1).*H(n,b.*k0).*J(n,a.*k2)-k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,b.*k2).*J(n,a.*k2)+k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,a.*k2).*J(n,b.*k2)-k1.*k2.*DH(n,b.*k0).*DH(n,a.*k2).*J(n,a.*k1).*J(n,b.*k2))))./J(n,a.*k1);
C21=-((k2.*(DJ(n,b.*k0).*H(n,b.*k0)-DH(n,b.*k0).*J(n,b.*k0)).*(k2.*DJ(n,a.*k1).*H(n,a.*k2)-k1.*DH(n,a.*k2).*J(n,a.*k1)))./(I.^n.*(-(k0.*k2.*DJ(n,a.*k1).*DJ(n,b.*k2).*H(n,b.*k0).*H(n,a.*k2))-k0.*k1.*DH(n,b.*k2).*DJ(n,a.*k2).*H(n,b.*k0).*J(n,a.*k1)+k0.*k1.*DH(n,a.*k2).*DJ(n,b.*k2).*H(n,b.*k0).*J(n,a.*k1)+k1.*k2.*DH(n,b.*k0).*DJ(n,a.*k2).*H(n,b.*k2).*J(n,a.*k1)+k0.*k2.*DH(n,b.*k2).*DJ(n,a.*k1).*H(n,b.*k0).*J(n,a.*k2)-k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,b.*k2).*J(n,a.*k2)+k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,a.*k2).*J(n,b.*k2)-k1.*k2.*DH(n,b.*k0).*DH(n,a.*k2).*J(n,a.*k1).*J(n,b.*k2))));
D21=(k2.*(DJ(n,b.*k0).*H(n,b.*k0)-DH(n,b.*k0).*J(n,b.*k0)).*(k2.*DJ(n,a.*k1).*H(n,a.*k2)-k1.*DH(n,a.*k2).*J(n,a.*k1)).*(k1.*DJ(n,a.*k2).*J(n,a.*k1)-k2.*DJ(n,a.*k1).*J(n,a.*k2)))./(I.^n.*(-(k2.*DJ(n,a.*k1).*H(n,a.*k2))+k1.*DH(n,a.*k2).*J(n,a.*k1)).*(-(k0.*k2.*DJ(n,a.*k1).*DJ(n,b.*k2).*H(n,b.*k0).*H(n,a.*k2))-k0.*k1.*DH(n,b.*k2).*DJ(n,a.*k2).*H(n,b.*k0).*J(n,a.*k1)+k0.*k1.*DH(n,a.*k2).*DJ(n,b.*k2).*H(n,b.*k0).*J(n,a.*k1)+k1.*k2.*DH(n,b.*k0).*DJ(n,a.*k2).*H(n,b.*k2).*J(n,a.*k1)+k0.*k2.*DH(n,b.*k2).*DJ(n,a.*k1).*H(n,b.*k0).*J(n,a.*k2)-k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,b.*k2).*J(n,a.*k2)+k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,a.*k2).*J(n,b.*k2)-k1.*k2.*DH(n,b.*k0).*DH(n,a.*k2).*J(n,a.*k1).*J(n,b.*k2)));
C01=(-J(n,b.*k0)+(k2.*H(n,b.*k2).*(DJ(n,b.*k0).*H(n,b.*k0)-DH(n,b.*k0).*J(n,b.*k0)).*(k2.*DJ(n,a.*k1).*H(n,a.*k2)-k1.*DH(n,a.*k2).*J(n,a.*k1)).*(k1.*DJ(n,a.*k2).*J(n,a.*k1)-k2.*DJ(n,a.*k1).*J(n,a.*k2)))./((-(k2.*DJ(n,a.*k1).*H(n,a.*k2))+k1.*DH(n,a.*k2).*J(n,a.*k1)).*(-(k0.*k2.*DJ(n,a.*k1).*DJ(n,b.*k2).*H(n,b.*k0).*H(n,a.*k2))-k0.*k1.*DH(n,b.*k2).*DJ(n,a.*k2).*H(n,b.*k0).*J(n,a.*k1)+k0.*k1.*DH(n,a.*k2).*DJ(n,b.*k2).*H(n,b.*k0).*J(n,a.*k1)+k1.*k2.*DH(n,b.*k0).*DJ(n,a.*k2).*H(n,b.*k2).*J(n,a.*k1)+k0.*k2.*DH(n,b.*k2).*DJ(n,a.*k1).*H(n,b.*k0).*J(n,a.*k2)-k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,b.*k2).*J(n,a.*k2)+k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,a.*k2).*J(n,b.*k2)-k1.*k2.*DH(n,b.*k0).*DH(n,a.*k2).*J(n,a.*k1).*J(n,b.*k2)))-...
    (k2.*(DJ(n,b.*k0).*H(n,b.*k0)-DH(n,b.*k0).*J(n,b.*k0)).*(k2.*DJ(n,a.*k1).*H(n,a.*k2)-k1.*DH(n,a.*k2).*J(n,a.*k1)).*J(n,b.*k2))./(-(k0.*k2.*DJ(n,a.*k1).*DJ(n,b.*k2).*H(n,b.*k0).*H(n,a.*k2))-k0.*k1.*DH(n,b.*k2).*DJ(n,a.*k2).*H(n,b.*k0).*J(n,a.*k1)+k0.*k1.*DH(n,a.*k2).*DJ(n,b.*k2).*H(n,b.*k0).*J(n,a.*k1)+k1.*k2.*DH(n,b.*k0).*DJ(n,a.*k2).*H(n,b.*k2).*J(n,a.*k1)+k0.*k2.*DH(n,b.*k2).*DJ(n,a.*k1).*H(n,b.*k0).*J(n,a.*k2)-k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,b.*k2).*J(n,a.*k2)+k2.^2.*DH(n,b.*k0).*DJ(n,a.*k1).*H(n,a.*k2).*J(n,b.*k2)-k1.*k2.*DH(n,b.*k0).*DH(n,a.*k2).*J(n,a.*k1).*J(n,b.*k2)))./(I.^n.*H(n,b.*k0));

function [out] = J(n, z)
out = besselj(n, z);

function [out] = H(n, z)
out = besselh(n, 2, z);

function [out] = DJ(n, z)
out = (1/2)*(J(-1+n,z)-J(1+n,z));

function [out] = DH(n, z)
out = (1/2)*(H(-1+n,z)-H(1+n,z));