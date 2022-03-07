clc
clear
N=input('Input the number of interval:');
h=1/N;
x=linspace(0,1,N+1);
f=inline('2','m');
u=inline('v*(1-v)','v');
for i=1:N+1
uex(i)=u(x(i));
end
for j=1:N
   Elem(j,:)=[j,j+1]; 
end
A=sparse(N+1,N+1);
L=sparse(N+1,1);
Db=[1,N+1];
for j=1:N
    A(Elem(j,:),Elem(j,:))=A(Elem(j,:),Elem(j,:))+[1/h -1/h;-1/h 1/h];
    L(Elem(j,:),1)= L(Elem(j,:),1)+h*f(x(Elem(j,1))/2+x(Elem(j,2))/2)*[1/2;1/2];
end
fullnodes=[1:N+1];
freenodes=setdiff(fullnodes,Db);
uh=zeros(N+1,1);
uh(freenodes,1)=A(freenodes,freenodes)\L(freenodes,1);
plot(x,uh,'r*',x,uex,'b');