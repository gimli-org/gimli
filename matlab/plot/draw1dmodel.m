function draw1dmodel(rho,thk,col,xlab,xlin,ylin,xl)

%% draw1dmodel(rho,thk,col,xlabel,xlin,xlim)

if nargin<5, ylin=0; end
if nargin<5, xlin=0; end
if nargin<4, xlab='\rho in \Omegam'; end
if nargin<3, col='r'; end

nl=length(rho);
if ylin, z0=0;
else z0=10.^floor(log10(thk(1)*0.9)); end
z=[z0;cumsum(thk(:))];z(end+1)=z(end)*2;
iz=1;ir=[];
for i=1:nl,
    ir(end+1:end+2)=i;
    iz(end+1:end+2)=i+1; 
end
iz(end)=[];
if ylin,
    if xlin, plot(rho(ir),z(iz),[col '-'],'LineWidth',1);
    else semilogx(rho(ir),z(iz),[col '-'],'LineWidth',1); end    
else
    if xlin, semilogy(rho(ir),z(iz),[col '-'],'LineWidth',1);
    else loglog(rho(ir),z(iz),[col '-'],'LineWidth',1); end
end
axis ij;grid on;
xlabel(xlab);ylabel('depth in m');
set(gca,'XAxisLocation','top');
if nargin>6, xlim(xl); end
xt=get(gca,'XTick');xtl=num2strcell(rndig(xt));
set(gca,'XTick',xt,'XTickLabel',xtl,'XTickMode','manual','XTickLabelMode','manual');
yt=get(gca,'YTick');ytl=num2strcell(rndig(yt));
set(gca,'YTick',yt,'YTickLabel',ytl,'YTickMode','manual','YTickLabelMode','manual');
