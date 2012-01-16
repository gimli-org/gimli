function [x, w] = gaulag(n)
% GAULAG - Gauss-Laguerre Integration Points
% [x,w] = gaulag(n)
% Given alf = 0.0, the parameter alpha of the Laguerre polynomials, this routine
% returns arrays x[1..n] and w[1..n] containing the abscissas and weights
% of the n-point Gauss-Laguerre quadrature formula. The smallest abscissa
% is returned in x[1], the largest in x[n].

% For a description of the following routines see
% Numerical Recipes, Press et a

EPS = 3.0e-11;
MAXIT = 10;

z = 0;
for i=1:n		%Loop over desired roots
	if i == 1
		z = 3.0/(1.0+2.4*n);
	elseif i == 2
			z = z + 15.0/(1.0+2.5*n);
	else
			ai = i-2;
			z = z + (1.0+2.55*ai)/(1.9*ai)*(z-x(ai));
   end
   for its = 1:MAXIT
      p1 = 1.0;
      p2 = 0.0;
      for j = 1:n
         p3 = p2;
         p2 = p1;
         p1 = ((2*j-1-z)*p2-(j-1)*p3)/j;
      end
		pp = n*(p1-p2)/z;
		z1 = z;
		z  = z1-p1/pp;
      if(abs(z-z1) <= EPS) break;
      end
   end
	x(i) = z;
	w(i) = -1/(pp*n*p2);
end