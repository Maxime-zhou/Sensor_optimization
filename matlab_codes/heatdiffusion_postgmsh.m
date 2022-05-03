
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing phase: writes data on file for post-processing wih GMSH
% Performs different actions according to the value of post_flag
% Writes history of temepratures. GRADIENT OUPUT NOT DEVELOPED!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if post_flag==0                                       % opens history file
 
 output_counter=0;   
 fmidhist=fopen('input/out.msh','w');                 
 fprintf(fmidhist,'%s\n%s\n%s\n','$MeshFormat','2.2 0 8','$EndMeshFormat');
 
elseif  post_flag==1,                         

 output_counter=output_counter+1;                     % displacement output
 fprintf(fmid,'%s\n','$NodeData','1',...
              '"U"','1',num2str(output_counter),'3',num2str(output_counter),'1');
 fprintf(fmid,'%d\n',analysis.NN);
 for n=1:analysis.NN
  fprintf(fmid,'%d %E\n',n,nodes(n).U);   
 end
 fprintf(fmid,'%s\n','$EndNodeData');
  
elseif  post_flag==2,                                 % final operations
    
 Sdim=0;   
 fmid=fopen('input/opti.geo','w');                    % options file
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
 fclose(fmid);

 dos(['gmsh.exe ', mshfile,' input/out.msh input/opti.geo']);  % runs gmsh with a dos command

end 
