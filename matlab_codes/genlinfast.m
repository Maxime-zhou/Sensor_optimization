
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% FILE: genlin_fast.m
% Linear analyses: scalar, elastic; 2D, axi, 3D 
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
La(2)=struct('tag','AX_therm_','ndof',1,'DG',2,'Sdim',3);  % thermal AX
La(3)=struct('tag','3D_therm_','ndof',1,'DG',3,'Sdim',3);  % thermal 3D
La(4)=struct('tag','2A_solid_','ndof',2,'DG',2,'Sdim',3);  % solid 2D
La(5)=struct('tag','2S_solid_','ndof',2,'DG',2,'Sdim',3);  % solid 2D
La(6)=struct('tag','AX_solid_','ndof',2,'DG',2,'Sdim',4);  % solid AX
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem dimensionality and tags for elemental subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Atag=La(analysis.type).tag;
ndof=La(analysis.type).ndof; 
DG=La(analysis.type).DG; 

readgmsh  % input file is read   

%pre_fract

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defines numbering of unknown nodal values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neq=0;
for e=1:analysis.NE,                      % global numbering of unknowns
 type=elements(e).type;
 Etag=Le(type).tag;
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
% preliminary allocation of workspace with Morse stocking
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
   nm(node)=nm(node)+Le(type).c;             % nm is incr. on connec. nodes
  end
 end 
end

startrow=zeros(analysis.neq,1);              % index of row start
ncoeffs=0;
for n=1:analysis.NN,                         % loop over nodes
 for dir=1:ndof                              % loop over directions
  dof=nodes(n).dof(dir);                     % DOF associated to node and idir
  if dof>0,                                  % only if unknown nodal value
   startrow(dof)=ncoeffs+1;
   ncoeffs=ncoeffs+ndof*nm(n);               % increments number of coefficients
  end
 end
end
Cm=zeros(ncoeffs,1);                         % allocates Cm
Jm=ones(ncoeffs,1);                          % allocates Jm vector
Im=ones(ncoeffs,1);                          % allocates Im vector
F=zeros(analysis.neq,1);                     % rhs vector
lenrow=zeros(analysis.neq,1);                % dim of row
count_coeff=0;

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
 Dne=ndof*ne;                                 % number of nodal values in one surface el.
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
 mat=elements(e).mat;
 Ke=eval([Etag Atag 'Ke(Xe,material(mat,:))']);       
 Le0=find(Ge>0);                           % local numbering of unknowns
 Ie=Ge(Le0);                               % global numbering 
 for i1=1:length(Ie)
  p=Le0(i1);                               % local number of unknown
  I=Ie(i1);                                % global number of unknown
  inRow=startrow(I);                       % initial position of row I
  dimRow=lenrow(I);                        % length of row I
  for i2=1:length(Ie)                      % loop over all columns in row I
   q=Le0(i2);                              % loc num of unknown in col i2
   J=Ie(i2);                               % glob num of unknown in col i2
   ipos=find(Jm(inRow:inRow+dimRow-1)==J); % verifies if coeff already exists
   if isempty(ipos),                       % if not creates a new coefficient
    dimRow=dimRow+1;
    lenrow(I)=dimRow;
    Im(inRow+dimRow-1)=I;
    Jm(inRow+dimRow-1)=J;
    Cm(inRow+dimRow-1)=Ke(p,q);
    count_coeff=count_coeff+1;
   else                                   % otherwise adds the contribution
    Cm(inRow+ipos-1)=Cm(inRow+ipos-1)+Ke(p,q);
   end
  end
 end
 LeD=find(Ge<0);                            % local numbering of prescr. b.c 
 F(Ie)=F(Ie)-Ke(Le0,LeD)*Ue(LeD);           % assemblage of rhs
end

K=sparse(Im,Jm,Cm);
clear Im Jm Cm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage phase: nodal loads due to surface tractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for e=1:analysis.NL,                      % for each element in loads
 type=loads(e).type;                      % type of the element
 Etag=Le(type).tag;
 ne=Le(type).ne;   
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
 F(Ie)=F(Ie)+Fe(Le0);                        % adds to global rhs
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------------------
disp('Solving system')
%-----------------------------------------------------------------------------------------

U=K\F;                                       % solution of linear system
clear K F                                    % clears memory

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
 ne=Le(type).ne;                          % number of nodes in element
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
 Sn(:,icomp)=Sn(:,icomp)./counter;  % naive average of stresses
end

%-----------------------------------------------------------------------------------------
disp('Launching GMSH')
disp('...........................')
%-----------------------------------------------------------------------------------------

%post_fract

genlin_postgmsh

