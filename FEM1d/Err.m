function [H1e,L2e]=Err(u,uh,x,N)

H1e=0; L2e=0;
h=1/N;

for i=1:N,
    p1=x(i); p2=x(i+1); mp=(p1+p2)/2;
    Duh=uh(i)*-1/h+uh(i+1)*1/h;
    
    H1e=H1e+h/6*(4*(uxe(mp)-Duh)^2+(uxe(x(i))-Duh)^2+(uxe(x(i+1))-Duh)^2);
   
    uhmp=uh(i)/2+uh(i+1)/2;
    L2e=L2e+h/6*((uh(i)-u(i))^2+(uh(i+1)-u(i+1))^2+4*(ue(mp)-uhmp)^2);
end
H1e=sqrt(H1e); L2e=sqrt(L2e);