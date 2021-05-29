[X,Z,B]=get_circleB;
circle_B_surf(X,Z,B)
function [X,Z,B]=get_circleB()
R=0.1;%Բ���뾶/m
N=100;%һ�����ϲ�������
c=3;
x=linspace(-c*R,c*R,N);%x���������
z=linspace(-c*R,c*R,N);%z���������
Bx=zeros(length(z),length(x));%��ʼ�� ����y���� ����x����
%By=zeros(N,N);
Bz=zeros(N,N);
r3=@(a,x,z)(R^2+x.^2+z.^2-2*R*x.*cos(a)).^(-3/2);
%����������x,z�������x,z��������ͻ
for i=1:length(z)
   for j=1:length(x)
       r3t=@(a)r3(a,x(j),z(i));%Bx,By,Bz������r3t,����������ټ����� ��aΪ����
       funx=@(a,z)R*z.*cos(a).*r3t(a);
       %funy=@(a,z)R*z.*sin(a).*r3t(a);
       funz=@(a,x)R*(R-x.*cos(a)).*r3t(a);
       Bx(i,j)=integral(@(a)funx(a,z(i)),0,2*pi); 
       %By(i,j)=integral(@(a)funy(a,z(i)),0,2*pi); 
       Bz(i,j)=integral(@(a)funz(a,x(j)),0,2*pi); 
       B(i,j)=sqrt(Bx(i,j)^2+Bz(i,j)^2);
   end
end
B=B*1e-7;%u0/4pi
[X,Z]=meshgrid(x,z);

end

function circle_B_surf(X,Z,B)
close all
R=0.1;%Բ���뾶/m
subplot(121)
surf(X,Z,B)
colorbar
shading interp
xlabel('X/m')
ylabel('Z/m')
zlabel('B')
title('xozƽ���ϴŸ�Ӧǿ�ȴ�С')
subplot(122)
surf(X,Z,log(B))
colorbar
xlabel('X/m')
ylabel('Z/m')
zlabel('log(B)')
title('xozƽ���ϴŸ�Ӧǿ�ȴ�С����ֵ')
end