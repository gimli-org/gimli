function res = rms(soll,ist,lolo)

% ARMS - Absolute Root Mean Square
% result = arms(a,b)

res = sqrt( mean( ( a - b ) .^2 ) );