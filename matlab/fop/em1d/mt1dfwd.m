function [rhoa,phi]=mt1dfwd(r,d,T)

% WAITALG - Wait algorithm
% = waitalg(resistivities,thicknesses,periods)

mu0 = 4*pi*10^-7;

fr= 1./T;
freq = fr*2*pi;
k1 = sqrt(i*freq*mu0/r(1));

g = 1;
nlay=min(length(r),length(d)+1);
for k = nlay-1:-1:1;    
    k1 = sqrt(i*freq*mu0/r(k));
    k2 = sqrt(i*freq*mu0/r(k+1));
    g = (g.*k2+k1.*tanh(k1*d(k)))./(k1+g.*k2.*tanh(k1*d(k)));    
end
z = i*freq./(k1.*g);

rhoa = mu0./freq.*abs(z).^2;
phi = angle(z);
zs=imag(z)./freq;
hind=phi>pi/4;
tind=phi<pi/4;
rhoh(hind)=2.*rhoa(hind).*(cos(phi(hind))).^2;
rhot(tind)=rhoa(tind)./(2.*(sin(phi(tind))).^2);
hh(hind) = sqrt(rhoa(hind)./freq(hind)/mu0).*...
    (sin(phi(hind))-cos(phi(hind)));
tt(tind) = sqrt(1./(rhoa(tind).*freq(tind)*mu0)).*...
    (cos(phi(tind))-sin(phi(tind)));

if nargout<1,
    fig=figure(1);
    subplot(2,1,1)
    loglog(T,rhoa,'-r'...
        ,T(hind),rhoh(hind),'og',...
        T(tind),rhot(tind),'*b');
    title({'1D by Wait and 2 layer model'});
    xlabel('T [s]');
    ylabel('\rho_a [\Omega m]');
    grid on;
    
    subplot(2,1,2)
    semilogx(T,phi*180/pi);
    set(gca,'YDir','Normal');
    xlabel('T [s]');
    ylabel('\Phi [°]');
    grid on;
end
%if nargin<2, rhoa=rhoa.*cos(phi)+rhoa.*sin(phi)*j; end
