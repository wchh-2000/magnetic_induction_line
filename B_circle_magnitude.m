[X,Z,B]=get_circleB;
circle_B_surf(X,Z,B)
function [X,Z,B]=get_circleB()
R=0.1;%圆环半径/m
N=100;%一个轴上采样点数
c=3;
x=linspace(-c*R,c*R,N);%x方向采样点
z=linspace(-c*R,c*R,N);%z方向采样点
Bx=zeros(length(z),length(x));%初始化 行数y长度 列数x长度
%By=zeros(N,N);
Bz=zeros(N,N);
r3=@(a,x,z)(R^2+x.^2+z.^2-2*R*x.*cos(a)).^(-3/2);
%匿名函数里x,z与外面的x,z变量不冲突
for i=1:length(z)
   for j=1:length(x)
       r3t=@(a)r3(a,x(j),z(i));%Bx,By,Bz都含有r3t,先算出来减少计算量 以a为变量
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
R=0.1;%圆环半径/m
subplot(121)
surf(X,Z,B)
colorbar
shading interp
xlabel('X/m')
ylabel('Z/m')
zlabel('B')
title('xoz平面上磁感应强度大小')
subplot(122)
surf(X,Z,log(B))
colorbar
xlabel('X/m')
ylabel('Z/m')
zlabel('log(B)')
title('xoz平面上磁感应强度大小对数值')
end