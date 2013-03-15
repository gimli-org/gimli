function x = cglscdpt(A,b,lam,C,D,P,dx,x0,maxiter,silent)

%% CGLSCDP - least squares normal equations solver
%% x = cglscdp(A,y,lam,C,D,P,dx)
%% solves implicitly the equation
%% ( (DS)'*DS + lam*C'C ) x = (DS)'*D*y (-lam*C'C dx)

if nargin<3, error('Too less input arguments!'); end
[n,m] = size(A); 
if nargin<3, lam=1; end
if nargin<4, C=1; end %speye(n); end
if nargin<5, D=1; end %ones(m,1); end
if nargin<6, P=1; end
if (nargin<7)||(isequal(dx,0)), dx=zeros(n,1); end
if (nargin<8)||(isequal(x0,0)), x0=zeros(n,1); end
if nargin<9, maxiter=1000; end
if nargin<10, silent=0; end

if min(size(D))==1, D=spdiags(D(:),0,length(D)); end
L=C;    
% Prepare for CG iteration.
%     PI = P\speye(size(P,1));
PI=P';
su=sum(P,1);
we=ones(size(su));
fi=find(su);we(fi)=1./su(fi);
PI=spdiags(we(:),0,length(we),length(we))*PI;
x = PI*x0;
z = D*(b - ((P*x)'*A)');%A*(P*x) % residuum of unregularized equation
p = PI*(A*(D*z));%(z'*D*A*P)';
acc=1e-8;
abbr = p'*p*acc; % goal for norm(r)^2
p=p-PI*(L*(x0+dx))*lam; % residuum of normal equation
r = p;
normr2 = r'*r;
% Iterate.
j=0;
if ~silent, wb=waitbar(0,'CGLS CDP'); else wb=1.11; end
fort=0;oldf=0;
t0=clock;
while(normr2>abbr)
  j=j+1;  
  q = D*((P*p)'*A)';  %;D*(A*(P*p));
  normr2old=normr2;
  Pp=P*p;
  alpha = normr2/(q'*q+Pp'*(L*Pp)*lam);
  x  = x + alpha*p;
  z  = z - alpha*q;
  r = PI*(A*(D*z))-PI*(L*(P*x+dx))*lam;%(z'*D*A*P)'-PI*(L*(P*x+dx))*lam;
  normr2 = r'*r;
  beta = normr2/normr2old;
  p = r + beta*p;
  fort=1+log10(normr2/abbr)/log10(acc);
  if ~silent&&(fort>oldf+0.05),
    waitbar(fort,wb);
    oldf=fort;
  end
end
if ishandle(wb), close(wb); end
x=P*x;
if ~silent, message(sprintf('Solved weighted normal equations in %.1fs, %d iterations,lam=%.1f',...
    etime(clock,t0),j,lam)); end