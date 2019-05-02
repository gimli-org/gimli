function savevtkmesh(Mesh,filename,ipvalue)

% SAVEVTKMESH - Save Mesh in VTK format
% savevtkmesh(Mesh,filename[,ipvalue])

% Mesh=[];
if nargin<2, 
    if nargin<1, error('Specify mesh'); end
    if nargin<2, filename='test.vtk'; end
end
if nargin<3, ipvalue=0; end
fid=fopen(filename,'w');
fprintf(fid,'# vtk DataFile Version 3.0\r\n');
fprintf(fid,'Created by dcmatlab\r\nASCII\r\nDATASET UNSTRUCTURED_GRID\r\n');
fprintf(fid,'POINTS %d double\r\n',Mesh.nnodes);
if size(Mesh.node,2)==3,
    fprintf(fid,'%g\t%g\t%g\r\n',Mesh.node');
else
    mm=Mesh.node;
%     mm(:,3)=Mesh.node(:,end);mm(:,2)=0;
    mm(:,3)=0;
    fprintf(fid,'%g\t%g\t%g\r\n',mm');
end
fprintf(fid,'CELLS %d %d\r\n',Mesh.ncells,sum(Mesh.cellnodes)+Mesh.ncells);
c = 0;
Mesh.cell = reshape(Mesh.cell',1,numel(Mesh.cell));
for i=1:Mesh.ncells
    fprintf(fid,'%d\t',Mesh.cellnodes(i),Mesh.cell(c+1:c+Mesh.cellnodes(i)-1)-1);
    fprintf(fid,'%d\r\n',Mesh.cell(c+Mesh.cellnodes(i))-1);
    c = c+Mesh.cellnodes(i);
end
fprintf(fid,'CELL_TYPES %d\r\n',Mesh.ncells);
cellmap = [1 3 5 10 10 10 10 12];
celltype = cellmap(Mesh.cellnodes);
fprintf(fid,'%d ',celltype);
fprintf(fid,'\r\n');
fprintf(fid,'CELL_DATA %d\r\n',Mesh.ncells);
fprintf(fid,'SCALARS Resistivity double 1\r\n');
fprintf(fid,'LOOKUP_TABLE default\r\n');
fprintf(fid,'%g ',Mesh.cellattr);
fprintf(fid,'\r\n');
if min(Mesh.cellattr)>0, % also put in log10(resistivity)
    fprintf(fid,'SCALARS Resistivity(log10) double 1\r\n');
    fprintf(fid,'LOOKUP_TABLE default\r\n');
    fprintf(fid,'%g ',log10(Mesh.cellattr));
    fprintf(fid,'\r\n');
end
if isfield(Mesh,'cellattr2log10')&&~isempty(Mesh.cellattr)&&(Mesh.cellattr(1))-Mesh.cellattr2(1)>1e-3, %write cellattr2
    fprintf(fid,'SCALARS Resistivity(log10) double 1\r\n');
    fprintf(fid,'LOOKUP_TABLE default\r\n');
    fprintf(fid,'%g ',Mesh.cellattr2);
    fprintf(fid,'\r\n');
end
if (nargin<=2)&&isfield(Mesh,'cellattr3')&&(length(Mesh.cellattr3)==Mesh.ncells),
    ipvalue=Mesh.cellattr3;
end
if length(ipvalue)==Mesh.ncells, 
    fprintf(fid,'SCALARS IP double 1\r\n');
    fprintf(fid,'LOOKUP_TABLE default\r\n');
    fprintf(fid,'%g ',ipvalue);
    fprintf(fid,'\r\n');    
end
if isfield(Mesh,'cellattr4')&&(length(Mesh.cellattr4)==Mesh.ncells),
    fprintf(fid,'SCALARS imag_res double 1\r\n');
    fprintf(fid,'LOOKUP_TABLE default\r\n');
    fprintf(fid,'%g ',Mesh.cellattr4);
    fprintf(fid,'\r\n');        
end
fclose(fid);