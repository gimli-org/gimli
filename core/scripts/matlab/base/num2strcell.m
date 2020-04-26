function ce=num2strcell(vec,ff)

% NUM2STRCELL - converts numerical array in string cell
% cells = num2strcell(array)

if (nargin<2)||(~ischar(ff)), ff='%g'; end
ce={};
for i=1:length(vec),
    ce{i}=strrep(sprintf(ff,vec(i)),'e+00','e+');
end