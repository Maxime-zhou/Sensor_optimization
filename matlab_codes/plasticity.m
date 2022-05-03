
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% FILE: plasticity.m
% Plastic analyses for linear hardening Von-Mises constitutive law
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

flag_post=0;                             % prepares output file
plasticity_postgmsh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of analysis type and elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

La(4)=struct('tag','2A_solid_','ndof',2,'DG',2,'Sdim',4); % only plane-strain!

Le(1)=struct('tag','B2_','ne',2,'ng',1,'c',1);            % limited choice of elements 
Le(2)=struct('tag','T3_','ne',3,'ng',1,'c',1); 
Le(8)=struct('tag','B3_','ne',3,'ng',2,'c',2); 
Le(9)=struct('tag','T6_','ne',6,'ng',3,'c',3); 

Atag=La(analysis.type).tag;
ndof=La(analysis.type).ndof;
DG=La(analysis.type).DG; 
Sdim=La(analysis.type).Sdim;

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
 elements(e).ng=Le(type).ng;                 % number of gauss point
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
analysis.ncoeffs = ncoeffs; 

KEP=spalloc(analysis.neq,analysis.neq,ncoeffs);% allocates sparse matrix
FuD=zeros(analysis.neq,1);                   % allocates rhs vector for tractions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage phase: system matrix and nodal forces due to imposed displ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------------------
disp(['Number of elements: ' num2str(analysis.NE)])
disp(['Number of nodes:    ',num2str(size(nodes,2))]);
disp(['Assembling initial elastic stiffness matrix']);
percold=0;
%-----------------------------------------------------------------------------------------

for e=1:analysis.NE,                        % stiffness matrix assemblage
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
 Le0=find(Ge>0);                            % local numbering of unknowns
 Ie=Ge(Le0);                                % global numbering 
 KEP(Ie,Ie)=KEP(Ie,Ie)+Ke(Le0,Le0);         % matrix assemblage
 LeD=find(Ge<0);                            % local numbering of prescr. b.c 
 FuD(Ie)=FuD(Ie)-Ke(Le0,LeD)*Ue(LeD);       % assemblage of rhs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage phase: nodal loads due to surface tractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fext=zeros(analysis.neq,1);               % allocates rhs vector for tractions
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
 Fext(Ie)=Fext(Ie)+Fe(Le0);                % adds to global rhs
end 

R=-Fext-FuD;                               % right hand side
U=-KEP\R;                                  % elastic solution for lambda=1
residref=sqrt(R'*R);                       % reference residuum 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:analysis.NN                       % loop on nodes
 nodes(n).DU(1:ndof)=LambdaT(2)* ...
                     nodes(n).U(1:ndof);  % saves prescribed displacement
 for dir=1:ndof
  dof=nodes(n).dof(dir);
  if dof>0 
   nodes(n).DU(dir)=LambdaT(2)*U(dof);    % saves elastic solution
  end 
 end
 nodes(n).U(1:ndof)=zeros(1:ndof);        % initialize displacement
end

for e=1:analysis.NE                       % loop on element
 ng=elements(e).ng ;                      % number of Gauss point
 elements(e).pg=zeros(ng,1) ;             % acc. plastic strain at Gauss pts 
 elements(e).Sg=zeros(ng,Sdim);           % stresses at Gauss pts
 elements(e).Dpg=zeros(ng,1) ;            % Inc. acc. plas. strain Gauss pts
 elements(e).Sghat=zeros(ng,Sdim);        % new stresses at Gauss pts
end
Fint=zeros(analysis.neq,1);               % Initialize internal forces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time sequence analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numstep=length(LambdaT)-1;                % Number of steps

%---------------------------------------------------------------------------
disp('Elastoplastic analysis')
disp(['Analysis phase. Total steps: ', num2str(numstep)])
%---------------------------------------------------------------------------

for nstep=0:numstep-1,                    % loop over all load steps
 iter=0;
 toll=1.d-4; 
 resid=1;

 if nstep>0
  for n=1:analysis.NN                     % initialize DU
   nodes(n).DU(:)=nodes(n).DU(:)* ...
      (LambdaT(nstep+2)-LambdaT(nstep+1))/(LambdaT(nstep+1)-LambdaT(nstep));
  end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton procedure for one single step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 while resid > toll*residref,    
  iter=iter+1;
  KEP(:,:)=0.d0;                          % sets Tangent matrix to zero
  Fint(:)=0;                              % sets Fint to zero
  
  for e=1:analysis.NE,                    % assemblage of KEP and Fint
   type=elements(e).type;
   Etag=Le(type).tag;
   ne=Le(type).ne;   
   Dne=ndof*ne;                           % number of nodal values
   mat=elements(e).mat;
   Xe=zeros(ne,DG);
   Ge=zeros(Dne,1);
   DUe=zeros(Dne,1);
   pos=1;
   for n=1:ne
    node=elements(e).nodes(n);   
    Xe(n,:)=nodes(node).coor;               
    Ge(pos:pos+ndof-1)=nodes(node).dof;
    DUe(pos:pos+ndof-1)=nodes(node).DU;
    pos=pos+ndof;
   end
   Sg=elements(e).Sg;                     % stress (Gauss pts)
   pg=elements(e).pg;                     % accumulated plastic strain (Gauss pts)
   [KEPe,Finte,elements(e).Sghat,elements(e).Dpg]=...
        eval([Etag Atag 'KTePlast'...
             '(Xe, material(mat,:),DUe,Sg,pg)']);
   Le0=find(Ge>0);                        % local numbering of unknowns
   Ie=Ge(Le0);                            % global numbering 
   Fint(Ie)=Fint(Ie)+Finte(Le0);          % assembling vector of internal forces
   KEP(Ie,Ie)=KEP(Ie,Ie)+KEPe(Le0,Le0);   % matrix assemblage
  end
  
  R=-Fext*LambdaT(nstep+2)-Fint;           % residuum vector
  resid=sqrt(R'*R);
  DDu=-KEP\R;                             % global prediction with consistent KEP   

%---------------------------------------------------------------------------
  disp(['Step: ', num2str(nstep), '; Lambda: ', num2str(LambdaT(nstep+2)), ...
        '; Iter: ', num2str(iter), '; Residual: ' num2str(resid)]);
%---------------------------------------------------------------------------

  for n=1:analysis.NN                     % updates displacement increment   
   for dir=1:ndof
    dof=nodes(n).dof(dir);
    if dof>0 
     nodes(n).DU(dir)=nodes(n).DU(dir)+DDu(dof);
    end 
   end
  end
 end                                      % end of Newton's iterations 

 for n=1:analysis.NN
  nodes(n).U=nodes(n).U+nodes(n).DU;      % updates displacement
 end
 for e=1:analysis.NE
  elements(e).pg=elements(e).pg...        % updates accumulated
                       + elements(e).Dpg; % plastic strains (Gauss pts)
  elements(e).Sg=elements(e).Sghat;       % updates stress (Gauss pts)
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Post-processing phase: computation of nodal stresses via extrapolation
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %--------------------------------------------------------------------------
 disp('Post processing data')
 %--------------------------------------------------------------------------

 DOUT=Sdim+1;                      
 Outn=zeros(analysis.NN,DOUT);            % initializes output matrix
 counter=zeros(analysis.NN,1);            % initializes counter list
 for e=1:analysis.NE,                     % loop over elements
  type=elements(e).type;
  Etag=Le(type).tag;
  connec=elements(e).nodes;               % gets element connectivty 
  ng=Le(type).ng;                         % number of Gauss pts in element 
  Outg(1:ng,1:Sdim)=elements(e).Sg;       % stress components
  Outg(1:ng,Sdim+1)=elements(e).pg(1:ng); % accumulated plastic strain 
  Outne=eval([Etag 'g2n(Outg)']);         % extrapolates at nodes
  Outn(connec,:)=Outn(connec,:)+Outne;    % adds to connec nodes
  counter(connec)=counter(connec)+1;
 end
 for icomp=1:DOUT
  Outn(:,icomp)=Outn(:,icomp)./counter;   % naive average of stresses
 end

 flag_post=1;                             % prepares output file
 plasticity_postgmsh
 
 clear Outn counter
 
end                                       % end loop over time steps

clear KEP Fext FuD Fint R DDu             % clears memory
 
%-----------------------------------------------------------------------------------------
disp('Launching GMSH')
disp('...........................')
%-----------------------------------------------------------------------------------------

flag_post=2;                             % prepares output file
plasticity_postgmsh
