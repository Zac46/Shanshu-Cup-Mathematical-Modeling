clear;
syms ini_theta;
ini_function=(10+25*9.8)*sin(ini_theta)*0.8+(10+25*9.8)*cos(ini_theta)*0.05-200*ini_theta;
initial_theta=vpasolve( ini_function, ini_theta,[0.0000000001,2]);%[0.0000000001,2]给出解范围，避免算出竖直的解
initial_theta=double(initial_theta);%类型转换，直接给的话是sym（符号类型）

syms  theta t;
syms  d_theta;%使用d_theta表示theta对时间一阶导
syms  dd_theta;%使用dd_theta表示theta对时间二阶导
syms m R_1 R_2 K l g;
theta_1=atan(R_1/l);      %右边的theta02的角
d=sqrt(R_1^2+(l)^2);   %两圆心间距
T=1/2*( 2/5*R_2^2*m + m*(d^2) )*( d_theta)^2 ;
%动能，1/2* I* d_theta^2
h=sqrt(R_1^2+l^2)*cos(theta+theta_1);
%球相对转轴高度
V=m*g*h+1/2*K*theta^2;
%势能

L=T-V;%拉格朗日量

func1=dt_diff( diff(L,d_theta) )-diff(L,theta);%==0

func1=subs(func1,[m,R_1,R_2,K,l,g],[25,0.05,0.1,200,0.8,9.8]);

f1=solve(func1,dd_theta);
%给出dd_theta的表达式

global s1;
s1=char(f1);
s1 = strrep(s1,'d_theta','y(2)') ;
s1 = strrep(s1,'theta','y(1)') ;
%传参用

tspan=[0 10];
y0= [initial_theta 0];
ops=odeset('RelTol',1e-10,'AbsTol',1e-15);
[time,result]=ode45(@ode,tspan,y0,ops);
%解微分方程，time为时间轴，result为结果，第一列为theta，第二列为d_theta


function dy=ode(x,y)
global s1;
dy=zeros(2,1);
dy(1)=y(2);
dy(2)=eval(s1);
end
%导入微分方程，通过全局变量s1

function fun=dt_diff(fun0)
syms  theta t;
syms  d_theta;
syms  dd_theta;
fun=diff(fun0,theta)*d_theta+diff(fun0,d_theta)*dd_theta+...
    diff(fun0,t);
end
%求t的全导数函数


