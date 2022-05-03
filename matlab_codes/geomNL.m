
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% FILE: geomNL.m
% Nonlinear elastic analyses (Saint-Venant Kirchhoff law): 2D, axi, 3D 
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessor phase: reads input from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name='Chap8_buckling.m';       % example of nonlinear elasticity (buckling)
eval(['run input/',name]);     % input file is read at once  

%-----------------------------------------------------------------------------------------
disp('...........................')
disp(['Reading input file: ' name])
%-----------------------------------------------------------------------------------------

% prepares output file
flag_post=0;
geomNL_postgmsh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of analysis type and elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

La(4)=struct('tag','2A_solid_','ndof',2,'DG',2,'Sdim',4);  % solid plane-strain

Le(1)=struct('tag','B2_','ne',2,'ng',1,'c',1); 
Le(2)=struct('tag','T3_','ne',3,'ng',1,'c',1); 
Le(3)=struct('tag','Q4_','ne',4,'ng',4,'c',2); 
Le(4)=struct('tag','P4_','ne',4,'ng',4,'c',1); 
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

%-----------------------------------------------------------------------------------------
disp(['Number of elements: ' num2str(analysis.NE)])
disp(['Number of nodes: ' num2str(analysis.NN)])
percold=0;
%-----------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage phase: nodal loads due to surface tractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F=zeros(analysis.neq,1);                  % allocates rhs vector for tractions
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
 Le0=find(Ge>0);                          % finds position of unknowns
 Ie=Ge(Le0);                              % and their number
 F(Ie)=F(Ie)+Fe(Le0);                     % adds to global rhs
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution phase: initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:analysis.NN                           % loop on nodes
 nodes(n).UD(1:ndof)=nodes(n).U(1:ndof);      % save prescribed displacement
 nodes(n).U(1:ndof)=zeros(1:ndof);            % Initialize displacement
end

Fint=zeros(analysis.neq,1);                   % Initialize internal forces
KT=spalloc(analysis.neq,analysis.neq,ncoeffs);% allocates sparse matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numstep = length(LambdaT)-1;                % number of steps
hist.Out=zeros(numstep,2);                  % for post proc issues

%---------------------------------------------------------------------------
disp('Nonlinear elastic analysis')
disp(['Analysis phase. Total steps: ', num2str(numstep)])
%---------------------------------------------------------------------------

for step=1:numstep,                        % loop over all load steps
 Dlambda=LambdaT(step+1)-LambdaT(step);    % inc. of loading factor
 Fext=F;
 for n=1:analysis.NN
  nodes(n).U=nodes(n).U+ ...               % initialisation of displ increment 
                  Dlambda*nodes(n).UD;     % using prescribed displ
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton procedure for one single step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 iter=0;                                   % resets iteration counter
 resid=1;                                  % fictitious initial value of res
 residref=1;                               % fictitious reference value of res
 toll=1.d-6;                               % tolerance on residuum
 while resid > toll*residref,
  iter=iter+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemblage of tangent matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  KT(:,:)=0.d0;                           % sets Tangent matrix to zero
  Fint(:)=0;                              % sets Fint to zero
  for e=1:analysis.NE,                    % loop over elements
   type=elements(e).type;
   Etag=Le(type).tag;
   ne=Le(type).ne;   
   Dne=ndof*ne;                           % number of nodal values in one surface el.
   mat=elements(e).mat;
   Xe=zeros(ne,DG);
   Ge=zeros(Dne,1);
   Ue=zeros(Dne,1);
   pos=1;
   for n=1:ne
    node=elements(e).nodes(n);   
    Xe(n,:)=nodes(node).coor;               
    Ge(pos:pos+ndof-1)=nodes(node).dof;
    Ue(pos:pos+ndof-1)=nodes(node).U;     % gets leement displacements
    pos=pos+ndof;
   end
   [KTe,Finte]=eval([Etag Atag 'KTeNL'...    % elemental tangent stiff matrix
            '(Xe,material(mat,:),Ue)']);     % and vector of internal forces 
   Le0=find(Ge>0);                           % local numbering of unknowns
   Ie=Ge(Le0);                               % global numbering 
   Fint(Ie)=Fint(Ie)+Finte(Le0);             % assembling internal forces
   KT(Ie,Ie)=KT(Ie,Ie)+KTe(Le0,Le0);         % matrix assemblage
  end
    
  R=-Fext-Fint;                           % residuum vector
  Du=-KT\R;                               % prediction with tangent matrix
  resid=sqrt(R'*R);
  if iter==1 
   residref=resid; 
  end
  for n=1:analysis.NN
   for dir=1:ndof
    dof=nodes(n).dof(dir);
    if dof>0 
     nodes(n).U(dir)=nodes(n).U(dir)+Du(dof);
    end 
   end
  end
  
%---------------------------------------------------------------------------
  disp(['Step: ', num2str(step), ' Iter: ', num2str(iter), ' Displ: ' num2str(nodes(3).U(1)) ' Residuum: ' num2str(resid)]);
%---------------------------------------------------------------------------
 
 end                                      % end of Newton's iterations 
 
 flag_post=1;                             % wirtes on files for postproc
 geomNL_postgmsh

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % force computation
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 work=0.d0;
 for e=1:analysis.NE,                     % loop over elements
  type=elements(e).type;
  Etag=Le(type).tag;
  ne=Le(type).ne;                         % number of nodes in element
  Dne=ndof*ne;                            % number of nodal values in one surface el.
  Xe=zeros(ne,DG);
  Ue=zeros(Dne,1);
  pos=1;
  for n=1:ne
   node=elements(e).nodes(n);   
   Xe(n,:)=nodes(node).coor;              % creates element
   Ue(pos:pos+ndof-1)=nodes(node).U;      % fills element nodal displ
   pos=pos+ndof;
  end
  mat=elements(e).mat;
  [KTe,Finte]=eval([Etag Atag 'KTeNL'... 
             '(Xe, material(mat,:), Ue)']);
  work=work+Ue'*Finte;                     
 end
 
 Ud=nodes(4).U(2);                        % for postprocessing
 hist.Out(step,:)=[nodes(6).U(1),work/Ud]; 
 
end                                       % end loop over time steps

clear KT Fext FDu Fint R Du               % clears memory

%-----------------------------------------------------------------------------------------
disp('Launching GMSH')
disp('...........................')
%-----------------------------------------------------------------------------------------

flag_post=2;                              % final post operations
geomNL_postgmsh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing: compares with buckling force from beam theory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w')                     
title('Buckling of doubly clamped beam')
xlabel('Mid-span deflection','Fontsize',12)
ylabel('Compressive load','Fontsize',12)
hold on
plot(abs(hist.Out(:,1)),abs(hist.Out(:,2)),'ok-')    % plots numerical solution 

H=1;
L=20;
J=1/12*H^3;                                  % inertia modulus
force=4*pi^2*J/(L^2);                        % analytcal buckling force
hista=hist.Out;               
hista(:,2)=force;               
plot(abs(hista(:,1)),abs(hista(:,2)),'k-.')  % plots analytical solution
legend('Mesh1',...
	   'Theory',...
	   'Location','SouthEast') 
