function x = tlscg(A,b,jmax,D,P,x0)

% CGLSCDP - truncated least squares normal equations solver
% x = tlscg(A,y,jmax,D,P,x0)
% solves implicitly jmax dg steps of the normal equation
% ( (DSP)'*DSP ) (P\x) = (DSP)'*D*y
% x0 - starting vector

if nargin<3, error('Too less input arguments!'); end
[m,n] = size(A); 
if nargin<3, jmax=20; end
if nargin<4, D=1; end %ones(m,1); end
if nargin<5, P=1; end
if (nargin<6)||(isequal(x0,0)), x0=zeros(n,1); end

if min(size(D))==1, D=spdiags(D(:),0,length(D)); end
% Prepare for CG iteration.
nn=size(P,2);
PI=P';
su=sum(P,1);
we=ones(size(su));
fi=find(su);we(fi)=1./su(fi);
PI=spdiags(we(:),0,length(we),length(we))*PI;
x = PI*x0;
z = D*(b - A*(P*x)); % residuum of unregularized equation
p = (z'*D*A*P)';
acc=1e-7;
abbr = p'*p*acc; % goal for norm(r)^2
r = p;
normr2 = r'*r;
% Iterate.
j=0;
wb=waitbar(0,'TLSCG');
fort=0;oldf=0;
t0=clock;
while(j<jmax)
  j=j+1;  
  q = D*(A*(P*p));
  normr2old=normr2;
  alpha = normr2/(q'*q);
  x  = x + alpha*p;
  z  = z - alpha*q;
  r = (z'*D*A*P)';
  normr2 = r'*r;
  beta = normr2/normr2old;
  p = r + beta*p;
  fort=j/jmax;
  if fort>oldf+0.05,
    waitbar(fort,wb);
    oldf=fort;
  end
end
close(wb);
x=P*x;
message(sprintf('Iterated %d steps of the normal equations in %.1fs',...
    jmax,etime(clock,t0)));