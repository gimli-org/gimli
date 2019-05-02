function [J,resp] = AmplitudeJacobian(G,mj)

% AMPLITUDEJACOBIAN - Calculate (real) amplitude jacobian from complex
%                     SNMR kernel function and given water content
% J = AmplitudeJacobian(K,m)
% J = AmplitudeJacobian(K) % for constant m=1
% [J,response] = AmplitudeJacobian(K) % for constant m=1

if nargin<2, mj=ones(size(G,2),1); end
if size(G,2)~=length(mj), error('Matrix and vector sizes mismatch!'); end
J=zeros(size(G));
dd = G*mj;
for m=1:size(G,2)
    J(:,m) = (real(G(:,m)).*real(dd) + imag(G(:,m)).*imag(dd))./abs(dd);
end
resp=abs(dd);