
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% File post_fract_T
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

a=1.;                                        % crack semilength
H=5.;                                        % plate halfwidth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K postprocessing phase from displacements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coor_K=0;
node_K=0;
for i=1:analysis.NN,
 coor=nodes(i).coor(:);   
 if (abs(coor(2))<1.d-6)& ...                % finds nodes on x_2=0 ... 
         (coor(1)<(a-1.d-5))&(coor(1)>coor_K) 
  node_K=i;       
  coor_K=coor(1);
 end 
end

delta_disp=2*nodes(node_K).U(2);
d=a-coor_K;
K_H=sqrt(pi*a)* ...                               % analytical solution for V->infty
        sqrt(sec(pi*a/(2*H)))
K_disp=material(1,1)/(8*(1-material(1,2)^2)) ...  % numerical K from displ.
          *sqrt(2*pi/d)*abs(delta_disp)
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K postprocessing phase from gtheta method: implemented only for T3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ri=.1;                                       % defines outer circle for theta=1 region
Re=.5;                                       % and inner circle for theta=0 region 
theta=zeros(analysis.NN,2);                 
for i=1:analysis.NN,
 coor=nodes(i).coor(:);   
 dist=sqrt((coor(1)-a)^2+coor(2)^2);         % distance from crack tip
 if dist<(Ri-1.d-5),
  theta(i,1)=1.;      
 elseif dist<(Re-1.d-5),
  theta(i,1)=1-(dist-Ri)/(Re-Ri);            % imposes linear theta field
 end 
end


GB=0.;
for e=1:analysis.NE,                      
 type=elements(e).type;
 ne=Le(type).ne;                        % number of nodes in element
 Etag=Le(type).tag;
 Dne=ndof*ne;                                % number of nodal values in one surface el.
 Ue=zeros(Dne,1);
 Ge=zeros(Dne,1);
 pos=1;
 for n=1:ne
  node=elements(e).nodes(n);   
  Xe(n,:)=nodes(node).coor;                  % creates element
  Ue(pos:pos+ndof-1)=nodes(node).U;
  Ge(pos:pos+ndof-1)=theta(node,:);
  pos=pos+ndof;
 end
 mat=elements(e).mat;
 GBe=eval([Etag Atag 'gtheta(Xe,material(mat,:),Ue,Ge)']);       
 GB=GB+GBe;             
end
K_GB=sqrt(material(1,1)/...                  % extrapolation of K_GB
             (1-material(1,2)^2)*2*GB)


