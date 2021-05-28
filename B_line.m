function B_line
R=cal_R()
theta = -pi:pi/10:pi;
x = R'*cos(theta);
y = R'*sin(theta);

rr = sqrt(x.^2 + y.^2);
% fx=-y;
% fy=x;
fx = -y./rr;%��λ�� ��λʸ��(fx,fy)
fy = x./rr;
scale = 0.2;
quiver(x, y, fx, fy,scale)
xlabel('x/m')
ylabel('y/m')
title('��ά������')
grid on
axis equal

function R=cal_R
%�õ��ȼ��仯�Ĵų���Ӧ��R(���㵽ֱ�߾���)
I=1;
r=0.01:0.001:0.1; %���㵽ֱ�߾���r��1cm ��1dm������1mm
B=2*1e-7*I./r;%B=u0*I/(2*pi.*r) u0=4*pi*1e-7
n=6;%��n���Ÿ���
Blabel=linspace(min(B),max(B),n);%n���ȼ���B,���Blabel
R=zeros(1,n);%��ʼ��
for i=1:length(Blabel)
   id=find_nearest(B,Blabel(i));%������B���ҵ���ӽ�Blabel(i)�ģ�����B���±�
   R(i)=r(id);%B���±���r���±�һһ��Ӧ
end

function id=find_nearest(B,obj)%������B���ҵ���ӽ�obj���ģ������±�
diff=abs(B-obj);
[m,id]=min(diff);%������С��
