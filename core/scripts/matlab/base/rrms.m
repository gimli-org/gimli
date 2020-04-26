function res = rrms(a,b,dummy)

% RRMS - Relative Root Mean Square (in %)
% result = rrms(a,b)

res = sqrt( mean( ( ( a - b ) ./ a ).^2 ) ) * 100;