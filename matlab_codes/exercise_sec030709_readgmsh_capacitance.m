
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reads GMSH file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mshfile=['input/' file];
fmid=fopen(mshfile,'r');            % opens mesh file

tline = fgets(fmid);                % reads a series of lines in preamble
tline = fgets(fmid); 
tline = fgets(fmid); 
tline = fgets(fmid); 

tline = fgets(fmid);                % gets string with number of nodes NN
NN=sscanf(tline,'%d');              % reads NN from string
position = ftell(fmid);             % memorizes position in mesh file
nodenum=zeros(NN,1);
for i=1:NN,                         % loop lines of section $Nodes
 tline = fgets(fmid); 
 h=sscanf(tline,'%d %f %f %f')';    % h contains node number and 3 coords
 nodenum(i)=h(1);                   % saves the node number in nodenum 
end
maxnodenum=max(nodenum);
fseek(fmid, position, 'bof');

nodes=repmat(struct('nodes',[],...  % allocates the list nodes
                    'dof',[],...
                    'U',[]),1,maxnodenum);   
for i=1:NN,                         % loop to read the nodes
 tline = fgets(fmid); 
 h=sscanf(tline,'%d %f %f %f')';    % reads node number and 3 coordinates
 nodes(h(1)).coor=h(2:1+DG);        % saves coordinates
 nodes(h(1)).dof=zeros(1,ndof);     % init to zero
 nodes(h(1)).U=zeros(1,ndof);       % init to zero
end
analysis.NN=NN;                     % save number of nodes

nsolid=numel(solid(:,1));           % number of sets assoc. to materials
ndbc=0; ndbcn=0; ntbc=0;
if exist('dbc')
 ndbc=numel(dbc(:,1));              % number of sets with imposed dbc
else    
 dbc=0;
end 
if exist('tbc') 
 ntbc=numel(tbc(:,1));              % number of sets with imposed tbc
else
 tbc=0;   
end   
if exist('dbcn')  
 ndbcn=numel(dbcn(:,1));            % number of nodes with imposed dbc
end 

tline = fgets(fmid); 
tline = fgets(fmid); 

NE=0;
NL=0;
tline = fgets(fmid);                % string with the number of elements
num=sscanf(tline,'%d');             % reads the global number of elements 
eltype=zeros(num,2);
position = ftell(fmid);             % memorizes position in mesh file
for i=1:num,             
 tline = fgets(fmid);         
 h=sscanf(tline,'%d')';
 physet=h(4);                       % physical set
 pos=find(solid(:,1)==physet);
 if ~isempty(pos) 
  NE=NE+1; 
  eltype(i,:)=[1 pos];
  continue
 end
 pos=find(tbc(:,1)==physet);
 NL=NL+length(pos);
 eltype(i,:)=[2 0];
 pos=find(dbc(:,1)==physet);
 for iL=1:length(pos)
  for n=6:length(h)
   nodes(h(n)).dof(dbc(pos(iL),2))=-1;
   nodes(h(n)).U(dbc(pos(iL),2))=dbc(pos(iL),3);
  end
 end      
 pos=find(constr(:)==physet);
 for iL=1:length(pos)
  for n=6:length(h)
   nodes(h(n)).dof(1)=1;
  end
 end      
end

fseek(fmid, position, 'bof');

elements=repmat(struct('type',0,...   % allocates list of elements 
                       'solid',0,...
                       'nodes',[]),1,NE);     
loads=repmat(struct('type',0,...      % allocates list of loads
                    'nodes',[],...
                    'dir',0,...
                    'val',0),1,NL);     

NE=0;
NL=0;
for i=1:num,             
 tline = fgets(fmid);         
 h=sscanf(tline,'%d')';
 physet=h(4);                         % physical set
 pos=eltype(i,2);
 switch eltype(i,1)
 case 1,    
  NE=NE+1;
  elements(NE).type=h(2);             % element type 
  elements(NE).mat=solid(pos,2);      % element material 
  elements(NE).nodes=h(6:length(h));  % connectivity of the element
 case 2,  
  pos=find(tbc(:,1)==physet);
  for iL=1:length(pos)
   NL=NL+1;
   loads(NL).type=h(2);            
   loads(NL).nodes=h(6:length(h));
   loads(NL).dir=tbc(pos(iL),2);
   loads(NL).val=tbc(pos(iL),3);
  end
 end 
end
analysis.NE=NE;
analysis.NL=NL;

for j=1:ndbcn   
 nodes(dbcn(j,1)).dof(dbcn(j,2))=-1;
 nodes(dbcn(j,1)).U(dbcn(j,2))=dbcn(j,3);
end    

fclose(fmid);

clear eltype nodenum

