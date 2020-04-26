function plotmymatrix(A,x,y,islog,cax)

% PLOTMYMATRIX - plot a matrix using imagesc and label
% plotmymatrix(A,x,y[,islog,caxis]);

if nargin<5, cax=[min(A(:)) max(A(:))]; end
if nargin<4, islog=0; end
if islog,
    imagesc(log(abs(A)));
    alpha(double(A>0));
    cax=log(cax);
else
    imagesc(A);
end
caxis(cax);
cb=colorbar;
if islog,
    ct=rndig(exp(get(cb,'YTick')),3);
    set(cb,'YTickLabel',num2strcell(ct));
end
xt=get(gca,'XTick');
xtl=num2strcell(rndig(x(xt),3));
yt=get(gca,'YTick');
ytl=num2strcell(rndig(y(yt),3));
set(gca,'XTickLabel',xtl,'YTickLabel',ytl);