function mesh2vtk(filename,Mesh,att,name,att2,name2,att3,name3,att4,name4);

% MESH2VTK - Save Mesh with property to VTK file
% mesh2vtk(filename,Mesh,property[,name,property2,name2,...])

if nargin<2, error('mesh2vtk(filename,Mesh,[property[,name]])'); end
if nargin<4, name='attribute'; end
if nargin<6, name2='attribute2'; end
if nargin<8, name3='attribute3'; end
if nargin<10, name4='attribute4'; end

if nargin<3, att=Mesh.cellattr; end
if isstruct(filename)&&isstr(Mesh),
   dummy=filename;filename=Mesh;Mesh=dummy; 
end
ctype=[1 5 10]; %point triangle tetrahedron
fid=fopen(filename,'w');
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'3d refraction\nASCII\nDATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %d double\n',Mesh.nnodes);
% ss='%f';for i=2:Mesh.dim, ss=[ss '\t%f']; end
ss='%f';for i=2:3, ss=[ss '\t%f']; end
node=Mesh.node;
if Mesh.dim==2, node(:,3)=Mesh.node(:,2);node(:,2)=0; end
fprintf(fid,[ss '\n'],node');
fprintf(fid,'CELLS %d %d\n',Mesh.ncells,Mesh.ncells*(Mesh.dim+2));
cells=zeros(Mesh.ncells,Mesh.dim+2);
cells(:,1)=Mesh.dim+1;
cells(:,2:Mesh.dim+2)=Mesh.cell-1;
ss='%d';for i=1:Mesh.dim+1, ss=[ss '\t%d']; end
fprintf(fid,[ss '\n'],cells');
fprintf(fid,'CELL_TYPES %d\n',Mesh.ncells);
fprintf(fid,'%d ',ones(Mesh.ncells,1)*ctype(Mesh.dim));
fprintf(fid,'\n');
fprintf(fid,'CELL_DATA %d\nSCALARS %s double 1\nLOOKUP_TABLE default\n',...
    Mesh.ncells,name);
fprintf(fid,'%f ',att);fprintf(fid,'\n');
if nargin>4, % 2 properties specified
    fprintf(fid,'SCALARS %s double 1\nLOOKUP_TABLE default\n',name2);
    fprintf(fid,'%f ',att2);fprintf(fid,'\n');
end
if nargin>6, % 3 properties specified
    fprintf(fid,'SCALARS %s double 1\nLOOKUP_TABLE default\n',name3);
    fprintf(fid,'%f ',att3);fprintf(fid,'\n');
end
if nargin>8, % 4 properties specified
    fprintf(fid,'SCALARS %s double 1\nLOOKUP_TABLE default\n',name4);
    fprintf(fid,'%f ',att4);fprintf(fid,'\n');
end
fclose(fid);