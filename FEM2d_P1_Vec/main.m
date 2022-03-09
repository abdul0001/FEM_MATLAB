%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM for  - div . (u_x,u_y) + u  = f    in  \Omega
%                   (u_x, u_y). n = u_N  on  \partial\Omega_N
%                               u = g    on  \partial\Omega_D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
L=1;% domain is (o,L)^2
N=8;% number of sub intervals in each sides of the square
[Coord,Elem,Nb,Db]=crisscross(L,N);
FullNodes=[1:size(Coord,1)];
FreeNodes=setdiff(FullNodes, unique(Db)); % DOF
nE = size(Elem,1);
 c1 = Coord(Elem(:,1),:);
 d21 = Coord(Elem(:,2),:)-c1;
 d31 = Coord(Elem(:,3),:) -c1;
 %*** Vector of element areas 4*|T|
 area4 = 2*(d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));
 %*** Assembly of stiffness matrix
 I = reshape(Elem(:,[1 2 3 1 2 3 1 2 3])',9*nE,1);
 J = reshape(Elem(:,[1 1 1 2 2 2 3 3 3])',9*nE,1);
 a = (sum(d21.*d31,2)./area4)';
 b = (sum(d31.*d31,2)./area4)';
 c = (sum(d21.*d21,2)./area4)';
% M=ones(nE,1)*[1/6 1/12 1/12 1/12 1/6 1/12 1/12 1/12 1/6];
% M=(M.*area4)'/4;
 A = [-2*a+b+c;a-b;a-c;a-b;b;-a;a-c;-a;c];
 A = sparse(I,J,A(:));
% M = sparse(I,J,M(:));
b=sparse(size(Coord,1),1); % global load vector
% Assembly of load vector b
for j=1:size(Elem,1)
    b(Elem(j,:))=b(Elem(j,:))+det([1 1 1; Coord(Elem(j,:),:)'])*f(sum(Coord(Elem(j,:),:))/3)/6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann Conditions
if (~isempty(Nb))
   for j=1:size(Nb,1)
       b(Nb(j,:))= b(Nb(j,:))+norm(Coord(Nb(j,1),:)-Coord(Nb(j,2),:))*...
                                     u_N(Coord(Nb(j,1),:),Coord(Nb(j,2),:))/2;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh=zeros(length(FullNodes),1);
% Dirichlet Conditions
if (~isempty(Db))
    Dbnodes=unique(Db);
    for j=1:size(Dbnodes,1)
       uh(Dbnodes(j),1)=ue(Coord(Dbnodes(j),:));
    end
end
b=b-A*uh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the linear system
uh(FreeNodes)=A(FreeNodes,FreeNodes)\b(FreeNodes);
% Exact solution at the nodes
u=u_nodes(Coord);
figure(1)
show(Coord,Elem,uh,u)





