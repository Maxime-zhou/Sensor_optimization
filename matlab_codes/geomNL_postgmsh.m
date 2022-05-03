
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing phase: writes data on file for post-processing wih GMSH
% Performs different actions according to the value of post_flag
% Writes history of displacements. STRESS OUTPUT NOT IMPLEMENTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_post==0                                % if first call

 fmidhist=fopen('input/out.msh','w');          % opens output file for history
 fprintf(fmidhist,'%s\n%s\n%s\n','$MeshFormat','2.2 0 8','$EndMeshFormat');

elseif flag_post==1,                           % if standard output

 ndof=La(analysis.type).ndof;                  % displacement output
 fprintf(fmidhist,'%s\n','$NodeData','1',... 
               ['"U"'],'1',num2str(step),'3',num2str(step),'3');
 fprintf(fmidhist,'%d\n',analysis.NN);
 for n=1:analysis.NN
  Uout=zeros(1,3);
  Uout(1:ndof)=nodes(n).U;
  fprintf(fmidhist,'%d %E %E %E \n',n,Uout);
 end
 fprintf(fmidhist,'%s\n','$EndNodeData');
  
elseif flag_post==2,                           % final operations
    
 fclose(fmidhist);                             % closes history file  
 
 fmid=fopen('input/opti.geo','w');             % write options file
 fprintf(fmid,'%s\n','General.BackgroundGradient=0;');
 fprintf(fmid,'%s\n','Mesh.SurfaceEdges=0;');
 fprintf(fmid,'%s\n','Mesh.SurfaceFaces=0;');
 fprintf(fmid,'%s\n','View[0].Visible = 1;');
 ndof=La(analysis.type).ndof;
 fprintf(fmid,'%s\n','View[0].VectorType = 5;'); 
 DisplFactor=1;   
 fprintf(fmid,'%s\n','View[0].DisplacementFactor = ',num2str(DisplFactor),';');
 fprintf(fmid,'%s\n','View[0].NbIso=21;');
 fprintf(fmid,'%s\n','View[0].IntervalsType=3;');
 for n=1:ndof-1
  fprintf(fmid,'%s%d%s\n','View[',n,'].Visible=0;');
  fprintf(fmid,'%s%d%s\n','View[',n,'].NbIso=21;');
  fprintf(fmid,'%s%d%s\n','View[',n,'].IntervalsType=3;');
 end
 fclose(fmid);
 
 dos(['gmsh.exe ', mshfile,' input/out.msh input/opti.geo']);  % runs gmsh with a dos command

end 

