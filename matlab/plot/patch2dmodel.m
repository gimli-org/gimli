function [cmin,cmax]=patch2dmodel(x,z,dM,MAL,N,Alpha)

% PATCH2DMODEL - Draw 2D (grid) Model
% patch2dmodel(x,z,Model[,OPT,N,ALPHA])
% patch2dmodel(Model[,OPT,N,ALPHA])
% x,z..model block positions
% dM..Model(length(x)-1,length(z)-1) of model values
% OPT..structure of possible fields
%      xdir - direction of  x axis (0/1)
%      cauto - automatic coloring (0/1)
%      cmin/cmax - minimum/maximum color values
%      cmap  - colormap (0-6 default,b2r,hot,gray,jet,cool)
%      log   - logarithmic coloring
%      style - Color boxes(0), filled contours(1), interpolated(2)
%      high  - plot z-axis with factor x
%Alpha=[];

stdmal=struct('cauto',1,'cmin',100,'cmax',500,'cmap',0,'xdir',0,'style',0,'high',1);
if isstruct(x),
   if nargin>3, Alpha=MAL; else Alpha=1; end
   if nargin>2, N=dM; else N=[]; end
   if nargin>1, MAL=z; else MAL=stdmal; end
   z=x.z;
   dM=x.M;
   x=x.x;
else
    if nargin<6, Alpha=1; end
    if nargin<5, N=[]; end
    if nargin<4, MAL=stdmal; end
end
if ~isfield(MAL,'cauto')&isfield(MAL,'cmin')&isfield(MAL,'cmax'), MAL.cauto=0; end
if min(Alpha(:))<0, Alpha=1; end
if isfield(MAL,'style')&&(MAL.style>0),
    draw2dmodel(x,z,dM,MAL,Alpha);
    return;
end
xz=zeros(size(x));
if isfield(N,'topo'),
    xz=interp1(N.topo(:,1),N.topo(:,2),x,'linear','extrap');
    z=-z;
end
if isempty(Alpha), Alpha=1; end
clog=(min(dM(:))>0);
if isfield(MAL,'xdir'), xdir=MAL.xdir; else xdir=0; end
if isfield(MAL,'cmax'), cmax=MAL.cmax; else cmax=1; end
if isfield(MAL,'cmin'), cmin=MAL.cmin; else cmin=0; end
if isfield(MAL,'cmap'), cmap=MAL.cmap; else cmap=0; end
if isfield(MAL,'clog'), clog=MAL.clog; end
if isfield(MAL,'log'), clog=MAL.log; end
if isfield(MAL,'cauto'), cauto=MAL.cauto; else cauto=1; end
if isfield(MAL,'canot'), canot=MAL.canot; else canot=0; end
if isfield(MAL,'cbar'), cbar=MAL.cbar; else cbar=1; end
if isfield(MAL,'elec'), elec=MAL.elec; else elec=0; end
if isfield(MAL,'high'), high=MAL.high; else high=1; end
if isfield(MAL,'style'), style=MAL.style; else style=0; end
if isfield(MAL,'cont'), cont=MAL.cont; else cont=15; end
if isfield(MAL,'alpha')&&(MAL.alpha==0), Alpha=0; end
if (cauto==0)&&(clog==1), cmin=log10(cmin);cmax=log10(cmax); end
cont=30;
if max(dM(:))==min(dM(:)), style=0; end
% if min(dM(:))*max(dM(:))<0, % symmetric blue-white-red
%     cmax=max(abs(dM(:)));
%     cmin=-cmax;
%     cmap=2;
%     clog=0;
%     cauto=0;
% end

if clog==1,
    dm=log10(reshape(dM,length(x)-1,length(z)-1))';
else
    dm=reshape(dM,length(x)-1,length(z)-1)';
end
uni=unique(dm(:));
if (length(uni)<10)&&(~isfield(MAL,'alfa')), Alpha=1; end
perc=5;
if isfield(MAL,'perc'), perc=MAL.perc; end
if cauto==1,
    if length(uni)<10,
        cmin=min(uni);
        cmax=max(uni);
    else
        [cmin,cmax]=interperc(dm(:),[perc 100-perc]);
%         [NN,VV]= hist(dm(:),100);
%         CN=cumsum(NN);
%         CN=CN/max(CN);
%         imin=max(find(CN<0.01));
%         imax=min(find(CN>0.99));
%         if isempty(imin), imin=1; end
%         if isempty(imax), imax=length(VV); end
%         cmin=VV(imin);
%         cmax=VV(imax);
        cmin=rndig(cmin);cmax=rndig(cmax);
    end
    if cmax<=cmin, cmax=max(dm(:));cmin=min(dm(:)); end
    if cmax<=cmin, cmin=cmax-0.01*abs(cmax); end
end
if isequal(size(Alpha),size(dM)), % Alpha shading
    if min(Alpha(:))<=0, Alpha(:)=1;MAL.noascal=1; end
    if isfield(MAL,'noascal')&&(MAL.noascal),
        AA=Alpha';
    else
        AA=log(Alpha');
        [nn,hh]=hist(AA(:),50);
        nnn=cumsum(nn)/numel(AA);
        mi=hh(min(find(nnn>0.1)));
        ma=hh(max(find(nnn<0.7)));
        AA=(AA-mi)/(ma-mi);
        AA(AA<0)=0;
        AA(AA>1)=1;
    end
else
    AA=ones(size(dM))';    
end
if style==1, %(>0) filled contours( or smooth surface)
    xx=(x(1:end-1)+x(2:end))/2;
    zz=(z(1:end-1)+z(2:end))/2;
    zz=[z(1);zz(:);z(end)];
    xx=[x(1);xx(:);x(end)];
    dm=[dm(1,:);dm;dm(end,:)];
    dm=[dm(:,1),dm,dm(:,end)];
    if style>1, %smooth surface
        pcolor(xx,zz,dm);
        shading flat
        grid off
        if size(Alpha)==size(dM),
            AA=[AA(1,:);AA;AA(end,:)];
            AA=[AA(:,1:end),AA(:,end),AA(:,end)];
            alpha(AA);
        end
        shading interp;
    else % filled contours
        contourf(xx,zz,dm,cont);
        shading flat
    end
else % colorboard
    switch cmap
        case 2, colmap=colormap(b2r);
        case 3, colmap=colormap(hot);
        case 4, colmap=colormap(gray);
        case 5, colmap=colormap(jet);
        case 6, colmap=colormap(cool);
        otherwise, colmap=colormap(jet);
    end
    if isfield(MAL,'creverse')&&(MAL.creverse>0), 
        colmap=flipud(colmap); 
        colormap(colmap);
    end
    lcm=length(colmap);
    ecolmap=colmap;
    if (length(x)<30)||(length(unique(dM))<8), ecolmap(:)=0.2; end
    if isfield(MAL,'showgrid'),
       if MAL.showgrid==1, ecolmap(:)=0.2; end
       if MAL.showgrid==2, ecolmap=colmap; end
    end
    cla reset;if cmin>=cmax, cmin=cmax*0.99-0.1; end
    caxis([cmin,cmax]);
    for i=1:size(dm,2),
        for k=1:size(dm,1),
            cind=getcindex(dm(k,i),cmin,cmax,lcm);
            col=colmap(cind,:)*AA(k,i)+1-AA(k,i);
            ecol=ecolmap(cind,:)*AA(k,i)+1-AA(k,i);
            patch(x([i i i+1 i+1]),z([k k+1 k+1 k])+xz([i i i+1 i+1]),col,'EdgeColor',ecol); 
%             patch(x([i i i+1 i+1]),z([k k+1 k+1 k])+xz([i i i+1 i+1]),...
%                 colmap(cind,:),'EdgeColor',ecolmap(cind,:),'FaceAlpha',AA(k,i),'EdgeAlpha',AA(k,i)); 
        end
    end
end
if high>0, 
%     axis equal; 
set(gca,'DataAspectRatio',[1 1 1]);
end
axis tight
set(gca,'XAxisLocation','top');
if ~isfield(N,'topo'), set(gca,'Ydir','reverse'); 
else set(gca,'Ydir','normal');  end
if xdir>0,
    set(gca,'XDir','reverse'); 
else
    set(gca,'XDir','normal');
end
set(gca,'XTickMode','auto','XTickLabelMode','auto');
set(gca,'YTickMode','auto','YTickLabelMode','auto');
set(gcf,'Renderer','zbuffer');
if cbar,
    hc=colorbar('horiz');%set(hc,'FontName','Symbol');
    set(hc,'DataAspectRatio',get(hc,'DataAspectRatio').*[1 32 1]);
    if(clog==1)
        clabel=rndig(10.^get(hc,'XTick'),3);
%         fi=find(clabel>10);
%         clabel(fi)=rndig(clabel(fi),3);
%         fi=find((clabel<=10)&(clabel>1));
%         clabel(fi)=round(10*clabel(fi))/10;
%         fi=find((clabel<1)&(clabel>=0.01));
%         clabel(fi)=round(100*clabel(fi))/100;
    else    
        clabel=rndig(get(hc,'XTick'),3);
    end
    xtl=num2strcell(clabel);
%     if ischar(canot), xtl{end-1}=canot; end
    %strrep(canot,'Ohm*m','WM'); end
    set(hc,'XTickMode','manual','XTickLabelMode','manual','XTickLabel',xtl);
end
if high>1,
    set(gca,'DataAspectRatio',[1 1/high 1]);
end
if length(z)>6,
    if size(z,1)==1,
        zz=[z(1) z(3) z(5:end)];
    else
        zz=[z(1);z(3);z(5:end)];
    end;
else
    zz=z;
end
zz=zz(:);
% xl=get(gca,'XTickLabel');
% if size(xl,1)>1, xl(end-1,1:3)='x/m'; end
% set(gca,'XTickLabel',xl);
xt=get(gca,'XTick');
xl=num2strcell(xt);
if length(xl)>2, xl{end-1}='x/m'; end
set(gca,'XTickLabel',xl,'XTickMode','manual','XTickLabelMode','manual');
% ddz=round(length(zz)/5);if ddz==0, ddz=1; end
% set(gca,'YTick',[zz(1:ddz:end-1);zz(end)]);
yt=get(gca,'YTick');
% yt(end:end+1)=[(yt(end-1)+yt(end))/2 yt(end)];
% set(gca,'YTick',yt);
% % yl=get(gca,'YTickLabel');yl(end-1,1:3)='z/m';yl(end-1,4:end)=' ';
yl=num2strcell(yt);
if isfield(N,'topo'), l=2; else l=length(yl)-1; end
yl{l}='z/m';set(gca,'YTickLabel',yl,'YTickMode','manual','YTickLabelMode','manual');
if elec,
    if isfield(N,'elec'),
        hold on
        if isfield(N,'topo')&&(~isempty(N.topo)),
            plot(N.elec(:,1),interp1(N.topo(:,1),N.topo(:,2),N.elec(:,1),'linear','extrap'),...
                'wv','MarkerEdgeColor','black','MarkerSize',5);
        else plot(N.elec(:,1),N.elec(:,2),'w.'); end
        hold off
    end
else
%     hold on
%     plot(x(1),0,'.','MarkerSize',1);
%     hold off
end
% if style==2, shading('interp'); end
% lx=length(x);lz=length(z);
box on;
% hold on;plot(x([1 lx lx 1 1]),z([1 1 lz lz 1]),'k-');hold off;
% set(line(x([1 lx lx 1 1]),z([1 1 lz lz 1])-0.01),'Color','black');
%Testversion n. 5 Zeilen einklammern
global libmmfile
if ~isequal(libmmfile,4),
    xl=xlim;zl=ylim;
    set(line(xl,zl),'Color','black');
    set(line(xl,fliplr(zl)),'Color','black');
    tv=[145 144 150 140 141 154 169 223 139 140 154 171];
    tt=text(mean(xl),mean(zl),char(255-fliplr(tv)));
    set(tt,'FontSize',24,'HorizontalAlignment','center','VerticalAlignment','middle');
end
if ischar(canot), 
    ax=gca;
    axes(hc);xli=get(gca,'Xlim');yli=get(gca,'Ylim');
    tt=text(xli(1),mean(yli),MAL.canot);
    set(tt,'VerticalAlignment','middle','HorizontalAlignment','right');
    axes(ax);
end %strrep(canot,'Ohm*m','WM'); end
if clog, cmin=10^cmin;cmax=10^cmax; end
if nargout<2, cmin=[cmin cmax]; end
