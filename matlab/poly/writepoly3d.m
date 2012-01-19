function writepoly3d(filename,Poly)

% WRITEPOLY3D - Write poly file for 3d mesher (tetgen)
% writepoly3d(filename,Poly)
% Poly...structure of 
%     .node (n by 4 array) x,y,z and marker (-99 electrodes)
%     .face cell array of index vectors into Poly.node
%     (note that .face must be 1-based while the file is 0-based)
%     .region (n by 4 array) x,y,z (optional marker and maximum cell size)
%     .hole (n by 4 array) x,y,z

% if size(Poly.node,2)<4, Poly.node(:,4)=0; end
nn=size(Poly.node,1);
fid=fopen(filename,'w');
fprintf(fid,'%d\t%d\t%d\t%d\n',nn,3,0,1);
node=ones(nn,5);
node(:,1)=(0:nn-1)';
node(:,2:size(Poly.node,2)+1)=Poly.node;
if isfield(Poly,'nodemarker'), node(:,5)=Poly.nodemarker; end
fprintf(fid,'%d\t%5f\t%5f\t%5f\t%d\n',node');
fprintf(fid,'%d\t%d\n',length(Poly.face),1);
facemarker=zeros(length(Poly.face));
if isfield(Poly,'facemarker')&&(length(Poly.facemarker)==length(facemarker)),
    facemarker=Poly.facemarker; end
for i=1:length(Poly.face),
   fprintf(fid,'1\t0\t%d\n',facemarker(i));
   fprintf(fid,'%d',length(Poly.face{i}));
   fprintf(fid,'\t%d',Poly.face{i}-1);
   fprintf(fid,'\n');
end
if isfield(Poly,'hole')&&(size(Poly.hole,1)>0),    
    fprintf(fid,'%d\n',size(Poly.hole,1)); %region marker
    hole=(0:size(Poly.hole,1)-1)';
    hole(:,2:4)=Poly.hole(:,1:3);
    fprintf(fid,'%d\t%10f\t%10f\t%10f\n',hole');
else
    fprintf(fid,'0\n'); %hole marker
end
if isfield(Poly,'region')&&(size(Poly.region,1)>0),    
    fprintf(fid,'%d\n',size(Poly.region,1)); %region marker
    region=(0:size(Poly.region,1)-1)';
    region(:,2:4)=Poly.region(:,1:3);
    region(1,6)=0;%default marker and attribute
    if size(Poly.region,2)>3, region(:,5)=Poly.region(:,4); end
    if size(Poly.region,2)>4, region(:,6)=Poly.region(:,5); end    
    fprintf(fid,'%d\t%10f\t%10f\t%10f\t%d\t%e\n',region');
else
    fprintf(fid,'0\n'); %region marker
end
fclose(fid);