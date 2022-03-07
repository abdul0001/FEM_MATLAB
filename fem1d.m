%%%%%%%%%%%%%%% FEM Code %%%%%%%%%%%%%%%%%%%%%%%
%%              -u"=f  in (0,1)
%%             u(0)=0; u(1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
for ip=1:1,
    
N=10*ip; % Mesh points
x=linspace(0,1,N+1);
h=1/N; di(ip)=1/N;
for i=1:N,
    Elem(i,:)=[i i+1];
end
Db=[1,N+1];
A=sparse(N+1,N+1);
L=sparse(N+1,1);
for i=1:N,
    A(Elem(i,:), Elem(i,:))=A(Elem(i,:),Elem(i,:))+[1/h -1/h;-1/h 1/h];
    L(Elem(i,:),1)=L(Elem(i,:),1)+h*f(x(Elem(i,1))/2+x(Elem(i,2))/2)*[1/2;1/2];
end
fullnodes=[1:(N+1)];
freenodes=setdiff(fullnodes,Db);
uh=zeros(N+1,1);
uh(freenodes)=A(freenodes,freenodes)\L(freenodes,1);
for j=1:(N+1),
    u(j,1)=ue(x(j));
end

[H1e(ip),L2e(ip)]=Err(u,uh,x,N);

end
for jp=1:(ip-1)
    ocH1(jp)=log(H1e(jp)/H1e(jp+1))/log(di(jp)/di(jp+1));
    ocL2(jp)=log(L2e(jp)/L2e(jp+1))/log(di(jp)/di(jp+1));
end

plot(x,uh,'b+',x,u,'r')


