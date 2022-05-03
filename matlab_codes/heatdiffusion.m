
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% FILE: heatdiffusion.m
% Linear analyses: 2D, axi, 3D 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessor phase: reads input from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[name,pname]=uigetfile('input/*.m',...   % gets name of input file
                       'Choose input file');

eval(['run input/',name]);               % input file is read at once  

%-----------------------------------------------------------------------------------------
disp('...........................')
disp(['Reading input file: ' name])
%-----------------------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of analysis type and elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

La(1)=struct('tag','2D_therm_','ndof',1,'DG',2,'Sdim',2);  % thermal 2D
La(2)=struct('tag','AX_therm_','ndof',1,'DG',2,'Sdim',2);  % thermal AXI
La(3)=struct('tag','3D_therm_','ndof',1,'DG',3,'Sdim',3);  % thermal 3D

Atag=La(analysis.type).tag;
ndof=La(analysis.type).ndof; 
DG=La(analysis.type).DG; 

Le(1)=struct('tag','B2_','ne',2,'c',1); 
Le(2)=struct('tag','T3_','ne',3,'c',1); 
Le(3)=struct('tag','Q4_','ne',4,'c',2); 
Le(4)=struct('tag','P4_','ne',4,'c',1); 
Le(8)=struct('tag','B3_','ne',3,'c',2); 
Le(9)=struct('tag','T6_','ne',6,'c',3); 

readgmsh  % input file is read   

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
M=spalloc(analysis.neq,analysis.neq,ncoeffs);% allocates sparse matrix
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
 perc=10*floor(10*e/analysis.NE);
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
 M(Ie,Ie)=M(Ie,Ie)+Me(Le0,Le0);             % matrix assemblage
 LeD=find(Ge<0);                            % local numbering of prescr. b.c 
 F(Ie)=F(Ie)-Ke(Le0,LeD)*Ue(LeD);           % assemblage of rhs
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
disp('Cholesky decomposition; implicit approach')
%-----------------------------------------------------------------------------------------

K=K+M/Dt;                                 
KC=chol(K);
clear K

%-----------------------------------------------------------------------------------------
disp('Time marching solution')
percold=0;
%-----------------------------------------------------------------------------------------

U=zeros(neq,1);
nstep=floor(tf/Dt);                          % number of time steps assuming constant dt
Dstep=floor(output_interval/Dt);             % every Dstep an output will be written on file
post_flag=0;
heatdiffusion_postgmsh                       % write initial conditions 

for step=1:nstep,                            % loop over time steps

%-----------------------------------------------------------------------------------------
 perc=10*floor(10*step/nstep);
 if perc~=percold,
  percold=perc;
  disp(['Time history: ' num2str(percold) '%'])
 end
%-----------------------------------------------------------------------------------------

 Fg=F+1/Dt*M*U;                              % implicit case
 U=KC\(KC'\Fg);

 if mod(step,Dstep)==0                       % every Dstep wirtes output
  for n=1:analysis.NN                        % loop over the nodes
   for dir=1:ndof                            % loop over all nodal values
    dof=nodes(n).dof(dir);                   % gets the global unknown num.
    if dof>0                                 % if it is unknwon
     nodes(n).U(dir)=U(dof);                 % fills nodal value from sol vector
    end 
   end
  end
  post_flag=1;
  heatdiffusion_postgmsh                     % writes temperatures to create movie
 end
 
end

clear KC M

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing phase: computation of nodal stresses via extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------------------
disp('Post processing data')
%-----------------------------------------------------------------------------------------

for n=1:analysis.NN             % loop over the nodes
 for dir=1:ndof                 % loop over all nodal values
  dof=nodes(n).dof(dir);        % gets the global unknown num.
  if dof>0                      % if it is unknwon
   nodes(n).U(dir)=U(dof);      % fills nodal value from sol vector
  end 
 end
end

Sdim=La(analysis.type).Sdim; 
Sn=zeros(analysis.NN,Sdim);             % initializes output matrix
counter=zeros(analysis.NN,1);           % initializes counter list
for e=1:analysis.NE,                    % computes nodal stresses

%-----------------------------------------------------------------------------------------
 perc=20*floor(5*e/analysis.NE);
 if perc~=percold,
  percold=perc;
  disp(['Post proc: ' num2str(percold) '%'])
 end
%-----------------------------------------------------------------------------------------

 type=elements(e).type;
 ne=Le(type).ne;                        % number of nodes in element
 Etag=Le(type).tag;
 Dne=ndof*ne;                           % number of nodal val in element
 Xe=zeros(ne,DG);                       % init matrix of coordinates
 Ue=zeros(Dne,1);                       % init list of nodal values
 connec=elements(e).nodes;              % gets element connectivty 
 pos=1;              
 for n=1:ne
  node=connec(n);   
  Xe(n,:)=nodes(node).coor;             % fills with coordinates of node n
  Ue(pos:pos+ndof-1)=nodes(node).U;     % fills nodal values of node n
  pos=pos+ndof;
 end
 mat=elements(e).mat;                   % gets element material
 Sg=eval([Etag Atag 'Sg(Xe,material(mat,:),Ue)']);  % at gauss points
 Sne=eval([Etag 'g2n(Sg)']);            % extrapolates at nodes
 Sn(connec,:)=Sn(connec,:)+Sne;         % adds to connec nodes
 counter(connec)=counter(connec)+1;
end   
for icomp=1:Sdim
 Sn(:,icomp)=Sn(:,icomp)./counter;      % naive average of stresses
end

%-----------------------------------------------------------------------------------------
disp('Launching GMSH')
disp('...........................')
%-----------------------------------------------------------------------------------------

post_flag=2;
heatdiffusion_postgmsh                  % post with movie of the temp field



