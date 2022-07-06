
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing phase: writes data on file for post-processing wih GMSH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmid=fopen('input/out.msh','w');                      % opens output file
fprintf(fmid,'%s\n%s\n%s\n','$MeshFormat','2.2 0 8','$EndMeshFormat');

% displacement output
for icomp=1:ndof                                      % for each component
 fprintf(fmid,'%s\n','$NodeData','1',...
              ['"U',num2str(icomp),'"'],'1','0.0','3','0','1');
 fprintf(fmid,'%d\n',analysis.NN);
 for n=1:analysis.NN                                  % for every node
  fprintf(fmid,'%d %E\n',n,nodes(n).U(icomp));   
 end
 fprintf(fmid,'%s\n','$EndNodeData');
end

% strain output
for icomp=1:Sdim
 fprintf(fmid,'%s\n','$NodeData','1',...
              ['"D',num2str(icomp),'"'],'1','0.0','3','0','1');
 fprintf(fmid,'%d\n',analysis.NN);
 for n=1:analysis.NN
  fprintf(fmid,'%d %E\n',n,Dn(n,icomp));   
 end
 fprintf(fmid,'%s\n','$EndNodeData');
end

% stress output
for icomp=1:Sdim
 fprintf(fmid,'%s\n','$NodeData','1',...
              ['"S',num2str(icomp),'"'],'1','0.0','3','0','1');
 fprintf(fmid,'%d\n',analysis.NN);
 for n=1:analysis.NN
  fprintf(fmid,'%d %E\n',n,Sn(n,icomp));   
 end
 fprintf(fmid,'%s\n','$EndNodeData');
end
fclose(fmid);

% writes options file
fmid=fopen('input/opti.geo','w');
fprintf(fmid,'%s\n','General.BackgroundGradient=0;');
fprintf(fmid,'%s\n','Mesh.SurfaceEdges=0;');
fprintf(fmid,'%s\n','Mesh.SurfaceFaces=0;');
fprintf(fmid,'%s\n','View[0].Visible = 1;');
fprintf(fmid,'%s\n','View[0].VectorType = 5;');
fprintf(fmid,'%s\n','View[0].NbIso=21;');
fprintf(fmid,'%s\n','View[0].IntervalsType=3;');
for n=1:ndof+Sdim-1
 fprintf(fmid,'%s%d%s\n','View[',n,'].Visible=0;');
 fprintf(fmid,'%s%d%s\n','View[',n,'].VectorType = 5;');
 fprintf(fmid,'%s%d%s\n','View[',n,'].NbIso=21;');
 fprintf(fmid,'%s%d%s\n','View[',n,'].IntervalsType=3;');
end

fclose(fmid);                                         % closes out.msh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% launches gmsh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dos(['gmsh.exe ', mshfile,' input/out.msh input/opti.geo']);  % runs gmsh with a dos command
% system(['gmsh', mshfile,' input/out.msh input/opti.geo']);  % runs gmsh with a dos command
