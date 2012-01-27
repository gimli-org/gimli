function [cmin,cmax]=tripatchmod(Mesh,att,alfa,MAL)%cmin,cmax)

% TRIPATCHMOD - Patch triangular model
% tripatchmod(Mesh) - shows only mesh
%             Mesh  .. triangular Mesh with nodes and cells
% tripatchmod(Mesh,field) - shows field on mesh
%                  field .. field to be used for color
% tripatchmod(Mesh,field,alpha) - shows field with alpha values
%                        alpha - transparency values
% tripatchmod(Mesh,field,options) or 
% tripatchmod(Mesh,field,alpha,options) - use options
%   options .. structure of patching options
%          .cauto - automatic colors
%          .cmin/.cmax - minimum/maximum colorbar
%          .clog - logarithmic colorbar
%          .cbar - draw colorbar (1=yes/0=no)
%          .cflip - flip colorbar upside-down
%          .canot - annotate colorbar (string)
%          .perc - percentage for automatic colorbar from 0 and 100 [5]
%          .xlim/ylim - limit of y axis

%  tripatchmod(Mesh,att,alfa);

cmin=0;
set(gca,'XTickMode','auto','XTickLabelMode','auto');
if nargin==0, error('No mesh specified!'); end
if nargin<4, MAL=[]; end
if ~isfield(MAL,'nocla')||(MAL.nocla==0), 
    cla reset; 
end
if (nargin>2)&&isstruct(alfa), MAL=alfa;alfa=ones(size(att)); end
if isfield(MAL,'clust')&&(MAL.clust>0),
%     MAL=struct('clog',0,'cbar',0,'cmap',7,'cauto',0,'cmin',0.5,'cmax',7.5,'clust',MAL.clust,'oldstyle',1);    
    MAL.clog=0;MAL.cbar=0;MAL.cmap=7;MAL.cauto=0;MAL.cmin=0.5;
    MAL.cmax=7.5;MAL.oldstyle=1;
end
if ~isfield(MAL,'oldstyle'), MAL.oldstyle=1; end
if (nargin<3)||(length(alfa)~=Mesh.ncells), alfa=ones(Mesh.ncells,1); end
if (nargin<2)||(length(att)~=Mesh.ncells), % pure mesh without colors
    patch('Vertices',Mesh.node,'Faces',Mesh.cell,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);
    axis equal tight
    if nargin>2,
        if isstruct(alfa), MAL=alfa; end
        if isfield(MAL,'xlim'), set(gca,'XLim',MAL.xlim); end
        if isfield(MAL,'ylim'), set(gca,'Ylim',MAL.ylim); end    
    end
    xtl=cellstr(get(gca,'XTickLabel'));
    xtl{end-1}='x/m';
    set(gca,'XTickMode','manual','XTickLabelMode','manual','XTickLabel',xtl);
    ytl=cellstr(get(gca,'YTickLabel'));
    ytl{end-1}='z/m';
    set(gca,'YTickMode','manual','YTickLabelMode','manual','YTickLabel',ytl);
    box on;
    set(gca,'XAxisLocation','top');    
    if isfield(Mesh,'boundmarker'),
        fb=find(Mesh.boundmarker>0);
        for i=1:length(fb),
            nn=Mesh.node(Mesh.bound(fb(i),:),:);
            line(nn(:,1),nn(:,2),'LineWidth',2,'Color','red');
        end
    end
    return;
    if isfield(Mesh,'cellattr'), att=Mesh.cellattr; 
    else att=(0:Mesh.ncells-1)'; end
%     else att=rand(Mesh.ncells,1)*2-1; end
end
nel=size(Mesh.cell,1);
if length(att)<nel, att(nel)=1; end
if length(att)>nel, att(nel+1:end)=[]; end
% old part
% if isfield(MAL,'cmin'), cmin=MAL.cmin; else cmin=min(att); end
% if isfield(MAL,'cmax'), cmax=MAL.cmax; else cmax=max(att); end
% if isfield(MAL,'clog'), islog=MAL.clog;
% else islog=(min(att)>0)&&(cmin>0); end
% if islog, 
%     att=log10(att); 
%     cmin=log10(cmin);
%     cmax=log10(cmax);
% end
% perc=5;if isfield(MAL,'perc'), perc=MAL.perc; end
% if (~isfield(MAL,'cauto'))||(MAL.cauto>0), 
%     if (length(unique(att))<10), cmin=min(att);cmax=max(att); 
%     else [cmin,cmax]=interperc(att,[perc 100-perc]); end
% end
% if (~isfield(MAL,'cmap')||(MAL.cmap~=0))&&((~islog)&&(max(att)*min(att)<0)),
%     cmap=colormap(b2r(64));
%     if ~isfield(MAL,'cauto')||(MAL.cauto==1), 
%         cmax=max(abs([cmin cmax]));cmin=-cmax; end
% else
%     if isfield(MAL,'cmap'),
%         if MAL.cmap==6, cmap=hsv(6); end
%         if MAL.cmap==7, cmap=[0 0 1;0 1 1;0 1 0;1 1 0;1 0 0;1 0 1;0 0 0]; end
%     else cmap=jet(64); end
% end
% if isfield(MAL,'cflip')&&(MAL.cflip), 
%     cmap=flipud(cmap);colormap(cmap);
% end
% if ~(cmax>cmin), cmax=cmin+1; end
% lcm=length(cmap)-1;
% emap=cmap;
% if (length(att)<100)||(length(unique(att))<2), emap(:)=0.2; end
% % emap(:)=0.2;%!!!
% for i=1:length(att),
%     cind=round(1+(att(i)-cmin)/(cmax-cmin)*lcm);
%     if cind<1, cind=1; end
%     if cind>lcm, cind=lcm; end
% %     patch(NODE(ELE(i,2:4),1),NODE(ELE(i,2:4),2),cmap(cind,:),'EdgeColor',emap(cind,:));
%     col=cmap(cind,:)*alfa(i)+1-alfa(i);
%     ecol=emap(cind,:)*alfa(i)+1-alfa(i);
%     set(patch(Mesh.node(Mesh.cell(i,1:Mesh.cellnodes(i)),1),Mesh.node(Mesh.cell(i,1:Mesh.cellnodes(i)),2),col,'EdgeColor',ecol),'LineStyle','none');
% end

[cols,cmin,cmax,cmap,islog]=getcols(att,MAL,alfa);
lstyle='none';
if unique(att)<10, lstyle='-'; end
if isfield(MAL,'oldstyle')&&(MAL.oldstyle>0),
    for i=1:size(cols,1), 
        cnodes=Mesh.cell(i,1:Mesh.cellnodes(i),:);        
        patch(Mesh.node(cnodes,1),Mesh.node(cnodes,2),cols(i,:),'LineStyle',lstyle,'EdgeColor','none'); 
%         patch(Mesh.node(Mesh.cell(i,1:3),1),Mesh.node(Mesh.cell(i,1:3),2),cols(i,:),'LineStyle','none'); 
    end
else
    p=patch('Vertices',Mesh.node,'Faces',Mesh.cell,'CData',cols,'FaceColor','flat','LineStyle','none');
end
axis equal tight
if isfield(MAL,'xlim'), set(gca,'XLim',MAL.xlim); end
if isfield(MAL,'ylim'), set(gca,'Ylim',MAL.ylim); end
caxis([cmin cmax]);
xtl=cellstr(get(gca,'XTickLabel'));
xtl{end-1}='x/m';
set(gca,'XTickMode','manual','XTickLabelMode','manual','XTickLabel',xtl);
ytl=cellstr(get(gca,'YTickLabel'));
ytl{end-1}='z/m';
set(gca,'YTickMode','manual','YTickLabelMode','manual','YTickLabel',ytl);
box on;
set(gca,'XAxisLocation','top');
if isfield(MAL,'XDir')&&isequal(MAL.XDir,'reverse'),
    set(gca,'XDir','reverse');
end
if isfield(MAL,'YDir')&&isequal(MAL.YDir,'reverse'),
    set(gca,'YDir','reverse');
end
if isfield(MAL,'plain')&&(MAL.plain>0),
    set(gca,'XTick',[],'YTick',[],'XColor','white','YColor','white');
    box off; end
if isfield(MAL,'xdir')&&(MAL.xdir>0),
   set(gca,'XDir','reverse'); end
%%
% clf;tripatchmod(Mesh);
if isfield(Mesh,'boundmarker'),
    fb=find(Mesh.boundmarker>0);
    for i=1:length(fb),
        nn=Mesh.node(Mesh.bound(fb(i),:),:);
        line(nn(:,1),nn(:,2),'LineWidth',1,'Color','black');
    end
end
%%
ax=gca;
if isfield(MAL,'FontSize'), set(ax,'FontSize',MAL.FontSize); end
if ~isfield(MAL,'cbar')||(MAL.cbar==1),
    colormap(cmap);
    cb=colorbar('horiz');%,'v6');
    dar=get(cb,'DataAspectRatio');
    set(cb,'DataAspectRatio',dar.*[1 32 1]);
    if islog,
        xt=get(cb,'XTick');
        xtl=num2strcell(rndig(10.^xt,3));
%         if isfield(MAL,'canot')&&ischar(MAL.canot), xtl{end-1}=MAL.canot; end
        set(cb,'XTickMode','manual','XTickLabelMode','manual','XTickLabel',xtl);
    end
    if isfield(MAL,'FontSize'), set(cb,'FontSize',MAL.FontSize); end
    if isfield(MAL,'canot'),
        %         set(cb,'YTick',mean(get(cb,'Ylim')),'YTickLabel',MAL.canot);
        axes(cb);xli=get(gca,'Xlim');yli=get(gca,'Ylim');
        tt=text(xli(1),mean(yli),MAL.canot);
        set(tt,'VerticalAlignment','middle','HorizontalAlignment','right');
        if isfield(MAL,'FontSize'), set(tt,'FontSize',MAL.FontSize); end
    end
end
if isfield(MAL,'cbar')&&(MAL.cbar==2),
    cb=colorbar;%,'v6');
    dar=get(cb,'DataAspectRatio');
    set(cb,'DataAspectRatio',dar.*[32 1 1]);
    if islog,
        yt=get(cb,'YTick');
        ytl=num2strcell(rndig(10.^yt,3));
%         if isfield(MAL,'canot')&&ischar(MAL.canot), xtl{end-1}=MAL.canot; end
        set(cb,'YTickLabelMode','manual','YTickLabel',ytl);
    end
    if isfield(MAL,'canot'),
        %         set(cb,'YTick',mean(get(cb,'Ylim')),'YTickLabel',MAL.canot);
        axes(cb);xli=get(gca,'Xlim');yli=get(gca,'Ylim');
        tt=text(mean(xli),yli(end),MAL.canot);
        set(tt,'VerticalAlignment','bottom','HorizontalAlignment','center');
    end
end
if islog, cmin=10^cmin;cmax=10^cmax; end
axes(ax);
