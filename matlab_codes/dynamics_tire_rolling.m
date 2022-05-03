
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Linear elastodynamics with central difference method: tire rolling down a slope
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessor phase: reads input from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name='Chap8_tire.m';                     
eval(['run input/',name]);                % input file is read at once  

%-----------------------------------------------------------------------------------------
disp('...........................')
disp(['Reading input file: ' name])
%-----------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of analysis type and elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

La(4)=struct('tag','2A_solid_','ndof',2,'DG',2,'Sdim',3);  % solid plane-strain
La(5)=struct('tag','2S_solid_','ndof',2,'DG',2,'Sdim',3);  % solid plane-stress
La(6)=struct('tag','AX_solid_','ndof',2,'DG',2,'Sdim',4);  % solid AXI
La(7)=struct('tag','3D_solid_','ndof',3,'DG',3,'Sdim',6);  % solid 3D

Atag=La(analysis.type).tag;
ndof=La(analysis.type).ndof; 
DG=La(analysis.type).DG; 

Le(1)=struct('tag','B2_','ne',2,'c',1); 
Le(2)=struct('tag','T3_','ne',3,'c',1); 
Le(3)=struct('tag','Q4_','ne',4,'c',2); 
Le(4)=struct('tag','P4_','ne',4,'c',1); 
Le(8)=struct('tag','B3_','ne',3,'c',2); 
Le(9)=struct('tag','T6_','ne',6,'c',3); 

readgmsh                                  % input file is read   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defines numbering of unknown nodal values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neq=0;
for e=1:analysis.NE,                      % global numbering of unknowns
 type=elements(e).type;
 ne=Le(type).ne;   
 connec=elements(e).nodes; 
 for n=1:ne
  node=connec(n);   
  for d=1:ndof,
   if nodes(node).dof(d)==0,              % if no dbc are enforced on node
    neq=neq+1;                            % increments equation number
    nodes(node).dof(d)=neq;               % and assigns value to current dof
   end  
  end
 end 
end
analysis.neq=neq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocation of sparse matrix (with estimate of the number of non zero coefficients)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nm=zeros(analysis.NN,1);                     % initialisation of nm
for e=1:analysis.NE,                         % loop over elements
 type=elements(e).type;
 connec=elements(e).nodes;                   % gets connectivity
 for n=1:Le(type).ne
  node=connec(n);   
  if nm(node)==0 
   nm(node)=nm(node)+Le(type).ne;
  else
   nm(node)=nm(node)+Le(type).c;             % nk is incr. on connec. nodes
  end
 end 
end
ncoeffs=ndof^2*sum(nm);                      % sum of all the terms in nm

K=spalloc(analysis.neq,analysis.neq,ncoeffs);% allocates sparse matrix
M=zeros(analysis.neq,1);                     % diagonal mass matrix
F=zeros(analysis.neq,1);                     % allocates rhs vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage phase: system matrix and nodal forces due to imposed displ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------------------
disp(['Number of elements: ' num2str(analysis.NE)])
percold=0;
%-----------------------------------------------------------------------------------------

for e=1:analysis.NE,                         % stiffness matrix assemblage
	
%-----------------------------------------------------------------------------------------
 perc=20*floor(5*e/analysis.NE);
 if perc~=percold,
  percold=perc;
  disp(['Assembling: ' num2str(percold) '%'])
 end
%-----------------------------------------------------------------------------------------

 type=elements(e).type;
 Etag=Le(type).tag;
 ne=Le(type).ne;   
 Dne=ndof*ne;                               % number of nodal values in one surface el.
 mat=elements(e).mat;
 Xe=zeros(ne,DG);
 Ge=zeros(Dne,1);
 Ue=zeros(Dne,1);
 pos=1;
 for n=1:ne
  node=elements(e).nodes(n);   
  Xe(n,:)=nodes(node).coor;                 % creates element
  Ge(pos:pos+ndof-1)=nodes(node).dof;
  Ue(pos:pos+ndof-1)=nodes(node).U;
  pos=pos+ndof;
 end
 Ke=eval([Etag Atag 'Ke(Xe,material(mat,:))']);       
 Me=eval([Etag Atag 'Me(Xe,material(mat,:))']);       
 Le0=find(Ge>0);                            % local numbering of unknowns
 Ie=Ge(Le0);                                % global numbering 
 K(Ie,Ie)=K(Ie,Ie)+Ke(Le0,Le0);             % matrix assemblage
 M(Ie)=M(Ie)+Me(Le0);                       % diagonal mass assemblage
 LeD=find(Ge<0);                            % local numbering of prescr. b.c 
 F(Ie)=F(Ie)-Ke(Le0,LeD)*Ue(LeD);           % assemblage of rhs

 Fe=eval([Etag Atag 'Fe(Xe,bf)']);          % body forces 
 F(Ie)=F(Ie)+Fe(Le0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage phase: nodal loads due to surface tractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for e=1:analysis.NL,                      % for each element in loads
 type=loads(e).type;                      % type of the element
 Etag=Le(type).tag;                       % element tag
 ne=Le(type).ne;                          % number of nodes 
 Dne=ndof*ne;                             % num of nodal values in one el
 idir=loads(e).dir;                       % traction direction
 val=loads(e).val;                        % traction value
 Xe=zeros(ne,DG);                        
 Ge=zeros(Dne,1);
 pos=1;
 for n=1:ne                                 
  node=loads(e).nodes(n);   
  Xe(n,:)=nodes(node).coor;               % matrix of nodal coords
  Ge(pos:pos+ndof-1)=nodes(node).dof;     % global numbering
  pos=pos+ndof;
 end
 Fe=eval([Etag Atag 'Fe_surf(Xe,val,idir)']); 
 Le0=find(Ge>0);                           % finds non-zero entries of dofe
 Ie=Ge(Le0); 
 F(Ie)=F(Ie)+Fe(Le0);                      % adds to global rhs
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time marching solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------------------
disp('Explicit approach: max stable dt')
%-----------------------------------------------------------------------------------------

uu=eigs(K,sparse(diag(M)),1);                % max stable time step (expl. analysis) 
stable_time_increment=2/sqrt(uu)

nstep=floor(tf/Dt);                          % number of time steps with constant dt
Dstep=floor(output_interval/Dt);             % every Dstep an output will be written on file

%-----------------------------------------------------------------------------------------
disp('Time marching solution')
percold=0;
%-----------------------------------------------------------------------------------------

U=zeros(neq,1);
Ud=zeros(neq,1);
Udd=M.\F;                                    % initial value of accelerations
Ud=Ud+Dt/2*Udd;                              % and of velocities

post_flag=0;
dynamics_tire_postgmsh                       % write initial conditions 

for step=1:nstep,                            % loop over time steps

%-----------------------------------------------------------------------------------------
 perc=10*floor(10*step/nstep);
 if perc~=percold,
  percold=perc;
  disp(['Time history: ' num2str(percold) '%  ' 'Displ: ' num2str(U(nodes(3).dof(2)))])
 end
%-----------------------------------------------------------------------------------------

% initial central difference operations
 U=U+Dt*Ud;                                  % actualisation of displacements
 Fg=F;
 % Fg=F-K*U;                                 % rhs without contact and small transf
 
% assemblage of internal forces
 for e=1:analysis.NE,                         
  type=elements(e).type;
  Etag=Le(type).tag;
  ne=Le(type).ne;   
  Dne=ndof*ne;                               % number of nodal values in one surface el.
  mat=elements(e).mat;
  Xe=zeros(ne,DG);
  Ge=zeros(Dne,1);
  Ue=zeros(Dne,1);
  pos=1;
  for n=1:ne
   node=elements(e).nodes(n);   
   Xe(n,:)=nodes(node).coor;                 % creates element
   Ge(pos:pos+ndof-1)=nodes(node).dof;
   pos=pos+ndof;
  end
  Ue=U(Ge);                                  % nodal values
  Fe=T3_2A_solid_Fint(Xe,material(mat,:),Ue);% internal nodal forces      
  Fg(Ge)=Fg(Ge)-Fe;                          % assemblage of rhs
 end
 
% explicit contact with penalty approach
 for n=1:analysis.NN                         % loop over the nodes
  dof=nodes(n).dof(1:2);
  if dof(1)>0
   nodes(n).U(1:2)=U(dof(1:2));
   x=nodes(n).coor(1:2)+U(dof(1:2))';
   if x(2)<x(1),                              
    force=penalty*(x(1)-x(2));
    Fg(dof(1))=Fg(dof(1))+(-1+f)*penalty*force;
    Fg(dof(2))=Fg(dof(2))+(1+f)*penalty*force;
   end       
  end 
 end

% resumes standard central difference scheme
 Udd=M.\Fg;                                  % acceleration
 Ud=Ud+Dt*Udd;                               % velocity
 
% writes output every Dstep steps
 if mod(step,Dstep)==0                       % writes on file displ for postprocessing  
  for n=1:analysis.NN                        % loop over the nodes
   for dir=1:ndof                            % loop over all nodal values
    dof=nodes(n).dof(dir);                   % gets the global unknown num.
    if dof>0                                 % if it is unknwon
     nodes(n).U(dir)=U(dof);                 % fills nodal value from sol vector
    end 
   end
  end
  post_flag=1;
  dynamics_tire_postgmsh                     % standard output
 end
  
end

clear KC M

post_flag=2;
dynamics_tire_postgmsh                       % final operations


