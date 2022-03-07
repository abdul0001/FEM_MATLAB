%%%%%%%%%%%%%%% FEM Code %%%%%%%%%%%%%%%%%%%%%%%
%%              -u"=f  in (0,1)
%%             u(0)=u(1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
N=100; % Mesh points
x=linspace(0,1,N+1);
h=1/N; 
for i=1:N
    Elem(i,:)=[i i+1];
end
Db=[1,N+1];
A=sparse(N+1,N+1);
L=sparse(N+1,1);
for i=1:N
    A(Elem(i,:), Elem(i,:))=A(Elem(i,:),Elem(i,:))+[1/h -1/h;-1/h 1/h];
    L(Elem(i,:),1)=L(Elem(i,:),1)+h*f(x(Elem(i,1))/2+x(Elem(i,2))/2)*[1/2;1/2];
end
fullnodes=[1:(N+1)];
freenodes=setdiff(fullnodes,Db);
uh=zeros(N+1,1);
uh(freenodes)=A(freenodes,freenodes)\L(freenodes,1);
for i=1:(N+1)
    u(i,1)=ue(x(i));
end
plot(x,uh,'b+',x,u,'r')


