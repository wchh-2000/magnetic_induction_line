function B_line
R=cal_R()
theta = -pi:pi/10:pi;
x = R'*cos(theta);
y = R'*sin(theta);

rr = sqrt(x.^2 + y.^2);
% fx=-y;
% fy=x;
fx = -y./rr;%单位化 单位矢量(fx,fy)
fy = x./rr;
scale = 0.2;
quiver(x, y, fx, fy,scale)
xlabel('x/m')
ylabel('y/m')
title('二维磁力线')
grid on
axis equal

function R=cal_R
%得到等间距变化的磁场对应的R(场点到直线距离)
I=1;
r=0.01:0.001:0.1; %场点到直线距离r从1cm 至1dm，步长1mm
B=2*1e-7*I./r;%B=u0*I/(2*pi.*r) u0=4*pi*1e-7
n=6;%画n条磁感线
Blabel=linspace(min(B),max(B),n);%n个等间距的B,组成Blabel
R=zeros(1,n);%初始化
for i=1:length(Blabel)
   id=find_nearest(B,Blabel(i));%从数组B中找到最接近Blabel(i)的，返回B内下标
   R(i)=r(id);%B的下标与r的下标一一对应
end

function id=find_nearest(B,obj)%从数组B中找到最接近obj数的，返回下标
diff=abs(B-obj);
[m,id]=min(diff);%差异最小的
