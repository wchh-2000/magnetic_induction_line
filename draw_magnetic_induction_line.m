[X,Z,A]=get_A();
circle_BA(X,Z,A)
function [X,Z,A]=get_A()
R=0.1;%圆环半径/m
N=500;%一个轴上采样点数
c=3;
x=linspace(-c*R,c*R,N);%x方向采样点
z=linspace(-c*R,c*R,N);%z方向采样点
[X,Z]=meshgrid(x,z);%矩阵X,Z
Xt=X;
X=abs(X);%A是关于|x|的函数
P=(X+R).^2+Z.*Z;%以下计算式中共有，先算出减少计算量 矩阵
M=4*R*X./P;%X,Z对应分量得到的m构成矩阵M, 椭圆模数k m=k^2
[K,E]=ellipke(M);%第一二类椭圆积分
A=sqrt(P)./X.*((X.^2+Z.^2+R^2)./P.*K-E);
X=Xt;%恢复X |x|用于算A X数据不应取绝对值
end

function circle_BA(X,Z,A)
close all
R=0.1;
F=abs(X).*A;

x0=R*[0,0.0947083,0.185857, 0.27102,0.350605,0.424542,0.493569, 0.557058,...
    0.613321,0.662027,0.70255, 0.737737,0.768563,0.79569,0.819633, 0.840806,...
    0.859552,0.87616,0.890875, 0.903911,0.915455,0.92567,0.934704, 0.942685,...
    0.949731,0.955943,0.961417, 0.966234,0.970469,0.97419,0.977455, 0.980317,0.982826,0.985021];
f0=funf0(x0);
figure(1);
hold on
contour(X,Z,F,f0);
plot([-R R],[0 0],'ro')
plot(-R,0,'r.')%垂直屏幕向外的电流
plot(R,0,'rx')%垂直屏幕向里的电流
grid on
xlabel('X/m')
ylabel('Z/m')
title('xoz平面上磁力线')

x1=R*[0.15,0.4,0.62,0.78,0.88,0.94];
f1=funf0(x1);
figure(2);
c=contour(X,Z,F,f1);%,'ShowText','on'
s=getcontourlines(c);%struct s .v值 .x,.y等值线上横纵坐标构成的数组
xloc=R*[0.4 0.8 0.8 0.9 1 1.02];%B的方向箭头的位置横坐标
for i=1:size(s,2)%s列数 为等值线条数 遍历等值线
   hold on
   xid=fix((i+1)/2);
   if s(i).x(1)>0%x正半轴范围
       xL=xloc(xid);%该等值线（磁力线）上B的方向箭头的位置横坐标
   else%x负半轴范围
       xL=-xloc(xid);
   end
   [m,id]=min(abs(s(i).x-xL));%离xL最近的点下标id
   vx=s(i).x(id); 
   vz=s(i).y(id);%B的方向箭头的位置坐标(vx,vz)
   [Bx,Bz,B]=cal_B(vx,vz);
   Dx=Bx/B;%方向向量单位化
   Dz=Bz/B;
   quiver(vx,vz,Dx,Dz,0.025,'r','MaxHeadSize',0.5)
end
plot([-R R],[0 0],'ro')
plot(-R,0,'r.')%垂直屏幕向外的电流
plot(R,0,'rx')%垂直屏幕向里的电流
set(gca,'XTick',[-0.3:0.02:0.3]);
xlabel('X/m')
ylabel('Z/m')
title('xoz平面上磁力线')
grid on 

color=['m' 'b' 'c' 'g' 'y' 'r'];
figure(3);
for i=1:size(s,2)%struct s 列数 为等值线条数
    hold on
    %旋转曲面
    n=25;%旋转一周所取点的个数
    theta = (0:n-1)/n*pi;
    X = s(i).x' * cos(theta);
    Y = s(i).x' * sin(theta);
    Z = s(i).y' * ones(1,n);%s(i).y对应z坐标
    ci=fix((i+1)/2);
    plot3(X,Y,Z,color(ci))
end
n=100;
t=linspace(0,2*pi,n);
x=R*cos(t);
y=R*sin(t);
z=zeros(1,n);
plot3(x,y,z,'k')
grid on
xlabel('X/m')
ylabel('Y/m')
zlabel('Z/m')
title('磁力线')
end

function f0=funf0(x0)%等值线标准值函数 x0为x轴上x=(-R,R)点集
R=0.1;
p=(x0+R).^2;
m=4*R*x0./(x0+R).^2;
[k,e]=ellipke(m);
f0=sqrt(p).*((x0.^2+R^2)./p.*k-e);
end

function [Bx,Bz,B]=cal_B(x,z)%由位置坐标x,z计算B
R=0.1;
r3=@(a)(R^2+x.^2+z.^2-2*R*x.*cos(a)).^(-3/2);
funx=@(a,z)R*z.*cos(a).*r3(a);
funz=@(a,x)R*(R-x.*cos(a)).*r3(a);
Bx=integral(@(a)funx(a,z),0,2*pi); 
Bz=integral(@(a)funz(a,x),0,2*pi);
B=sqrt(Bx^2+Bz^2);
end

function s = getcontourlines(c)
%  It takes the output of the contour function, and returns a struct array as output.
%  Each struct in the array represents one contour line. The struct has fields
% v, the value of the contour line   s(1).v
% x, the x coordinates of the points on the contour line s(1).x
% y, the y coordinates of the points on the contour line
sz = size(c,2); % Size of the contour matrix c
ii = 1; % Index to keep track of current location
jj = 1; % Counter to keep track of contour lines
while ii < sz % While we haven't exhausted the array
n = c(2,ii); % How many points in this contour?
s(jj).v = c(1,ii); % Value of the contour
s(jj).x = c(1,ii+1:ii+n); % X coordinates
s(jj).y = c(2,ii+1:ii+n); % Y coordinates
ii = ii + n + 1; % Skip ahead to next contour line
jj = jj + 1; % Increment number of contours
end
end