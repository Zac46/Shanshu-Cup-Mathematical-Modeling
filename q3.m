clear;
syms r_1 varphi_1 t;
syms dr_1 dvarphi_1;%dr_1和dvarphi_1分别为r_1和varphi_1的对时间的导数
syms ddr_1 ddvarphi_1;%ddr_1和ddvarphi_1分别为r_1和varphi_1的对时间的二阶导数
syms m1 m2 R_1 R_2 K l_1 l_2 omega f;
theta_1=atan(R_1/l_1);theta_2=atan(R_1/l_2);
r_2=5*r_1/6;
varphi_2=varphi_1+pi;
d_1=sqrt(R_1^2+l_1^2);d_2=sqrt(R_1^2+l_2^2);


T=m1/2*dr_1^2+m2/2*dt_diff(r_2)^2 ...
    + 1/2*( 2/5*m1*R_2^2+2/5*m2*R_2^2+m1*r_1^2 + m2*r_2^2)*dvarphi_1^2;
%动能项
theta= acos( ( d_1^2+d_2^2-( r_1+r_2 )^2 )/2/d_1/d_2 )-theta_2-theta_1 ;

V=K/2*(theta)^2;
%势能项



% V=K/2*( acos( ( d_1^2+d_2^2-( 11*r_1/6 )^2 )/2/d_1/d_2 )-theta_2-theta_1 )^2-...
%     r_1*cos(varphi_1)*5/11*f*cos(omega*t)-r_2*cos(varphi_1)*5/11*f*cos(omega*t);
%另外一种算法中的势能项，即将惯性力也看成一个含时的保守立场，解释较麻烦未采用



L=T-V;

func1=dt_diff( diff(L,dvarphi_1) )-diff(L,varphi_1)+5/6*f*cos(omega*t)*sin(varphi_1)*r_1;
func2=dt_diff( diff(L,dr_1) )-diff(L,r_1)-5/6*f*cos(omega * t)*cos(varphi_1);
%Euler-Lagrange方程

% func1=dt_diff( diff(L,dvarphi_1) )-diff(L,varphi_1);
% func2=dt_diff( diff(L,dr_1) )-diff(L,r_1);
%另外一种算法中的Euler-Lagrange方程，即将惯性力也看成一个含时的保守立场，解释较麻烦未采用


func1=subs(func1,[m1 m2 R_1 R_2 K l_1 l_2 omega f],[25 30 0.05 0.1 200 0.6 0.8 1.4 50]);
func2=subs(func2,[m1 m2 R_1 R_2 K l_1 l_2 omega f],[25 30 0.05 0.1 200 0.6 0.8 1.4 50]);
%带入相关常量


global s1 s2;
[f1,f2]=solve([func1,func2],[ddr_1 ddvarphi_1]);
s1=char(f1);
s2=char(f2);
s1 = strrep(s1,'dvarphi_1','y(4)') ;
s2 = strrep(s2,'dvarphi_1','y(4)') ;
s1 = strrep(s1,'varphi_1','y(2)') ;
s2 = strrep(s2,'varphi_1','y(2)') ;
s1 = strrep(s1,'dr_1','y(3)') ;
s2 = strrep(s2,'dr_1','y(3)') ;
s1 = strrep(s1,'r_1','y(1)') ;
s2 = strrep(s2,'r_1','y(1)') ;
%解出二阶导数满足关系

time=[];result=zeros(0,4);

tspan=[0 10];
y0=[sqrt(5)/10/11*6 -atan(2) 0 0];
% ops=[];
ops=odeset('RelTol',1e-11,'AbsTol',1e-15);
[t,yy]=ode45(@ode,tspan,y0,ops);
%解微分方程


end_time_index=length(yy(:,1));
while min(yy(1:length(yy(:,1)),1))-sqrt(5)/10/11*6<0
    for ind=1:length(yy(:,1))
        if yy(ind,1)-sqrt(5)/10/11*6<0
            end_time_index=ind;
            break;
        end
    end
        time=[time;t(1:end_time_index)];
        result=[result;yy(1:end_time_index,:)];
        y0=yy(end_time_index,:);y0(1)=sqrt(5)/10/11*6;
        y0(3)=-y0(3);tstart=t(end_time_index);
        tspan=[tstart 10];
        [t,yy]=ode45(@ode,tspan,y0,ops);
        
end
%计算碰撞相关代码，每次碰撞后以新状态为初始条件解微分方程

time=[time;t];
result=[result;yy];
%存放结果，result中从左到右分别为r_1,\varphi_1,dr_1,d\varphi_1




alpha = acos((d_2.^2 + (r_1+r_2).^2 -d_1.^2) / 2/ d_2 / (r_1+r_2) );
x=r_2 * cos(varphi_2) + d_2 * cos(varphi_1-alpha);
y=r_2 * sin(varphi_2) + d_2 * sin(varphi_1-alpha);
x_1=r_1*cos(varphi_1)-x;    y_1=r_1*sin(varphi_1)-y;
x_2=r_2*cos(varphi_2)-x;    y_2=r_2*sin(varphi_2)-y;
omega_1 = dt_diff( atan(y_1/x_1) );
omega_2 = dt_diff( atan(y_2/x_2) );



x_1=subs(x_1,[R_1 R_2 l_1 l_2],[0.05 0.1 0.6 0.8]);
y_1=subs(y_1,[R_1 R_2 l_1 l_2],[0.05 0.1 0.6 0.8]);
x_2=subs(x_2,[R_1 R_2 l_1 l_2],[0.05 0.1 0.6 0.8]);
y_2=subs(y_2,[R_1 R_2 l_1 l_2],[0.05 0.1 0.6 0.8]);
omega_1=subs(omega_1,[R_1 R_2 l_1 l_2],[0.05 0.1 0.6 0.8]);
omega_2=subs(omega_2,[R_1 R_2 l_1 l_2],[0.05 0.1 0.6 0.8]);
theta=subs(theta,[R_1 R_2 l_1 l_2],[0.05 0.1 0.6 0.8]);


string_x_1=char(x_1);   string_y_1=char(y_1);   
string_x_2=char(x_2);   string_y_2=char(y_2);
string_omega_1=char(omega_1);   string_omega_2=char(omega_2);
string_theta=char(theta);   

string_x_1 = strrep(string_x_1,'dvarphi_1','result(:,4)') ;
string_x_1 = strrep(string_x_1,'varphi_1','result(:,2)') ;
string_x_1 = strrep(string_x_1,'dr_1','result(:,3)') ;
string_x_1 = strrep(string_x_1,'r_1','result(:,1)') ;
string_x_1 = strrep(string_x_1,'^','.^') ;
string_x_1 = strrep(string_x_1,'*','.*') ;
string_x_1 = strrep(string_x_1,'/','./') ;

string_y_1 = strrep(string_y_1,'dvarphi_1','result(:,4)') ;
string_y_1 = strrep(string_y_1,'varphi_1','result(:,2)') ;
string_y_1 = strrep(string_y_1,'dr_1','result(:,3)') ;
string_y_1 = strrep(string_y_1,'r_1','result(:,1)') ;
string_y_1 = strrep(string_y_1,'^','.^') ;
string_y_1 = strrep(string_y_1,'*','.*') ;
string_y_1 = strrep(string_y_1,'/','./') ;

string_x_2 = strrep(string_x_2,'dvarphi_1','result(:,4)') ;
string_x_2 = strrep(string_x_2,'varphi_1','result(:,2)') ;
string_x_2 = strrep(string_x_2,'dr_1','result(:,3)') ;
string_x_2 = strrep(string_x_2,'r_1','result(:,1)') ;
string_x_2 = strrep(string_x_2,'^','.^') ;
string_x_2 = strrep(string_x_2,'*','.*') ;
string_x_2 = strrep(string_x_2,'/','./') ;

string_y_2 = strrep(string_y_2,'dvarphi_1','result(:,4)') ;
string_y_2 = strrep(string_y_2,'varphi_1','result(:,2)') ;
string_y_2 = strrep(string_y_2,'dr_1','result(:,3)') ;
string_y_2 = strrep(string_y_2,'r_1','result(:,1)') ;
string_y_2 = strrep(string_y_2,'^','.^') ;
string_y_2 = strrep(string_y_2,'*','.*') ;
string_y_2 = strrep(string_y_2,'/','./') ;

string_theta = strrep(string_theta,'dvarphi_1','result(:,4)') ;
string_theta = strrep(string_theta,'varphi_1','result(:,2)') ;
string_theta = strrep(string_theta,'dr_1','result(:,3)') ;
string_theta = strrep(string_theta,'r_1','result(:,1)') ;
string_theta = strrep(string_theta,'^','.^') ;
string_theta = strrep(string_theta,'*','.*') ;
string_theta = strrep(string_theta,'/','./') ;

string_omega_1 = strrep(string_omega_1,'dvarphi_1','result(:,4)') ;
string_omega_1 = strrep(string_omega_1,'varphi_1','result(:,2)') ;
string_omega_1 = strrep(string_omega_1,'dr_1','result(:,3)') ;
string_omega_1 = strrep(string_omega_1,'r_1','result(:,1)') ;
string_omega_1 = strrep(string_omega_1,'^','.^') ;
string_omega_1 = strrep(string_omega_1,'*','.*') ;
string_omega_1 = strrep(string_omega_1,'/','./') ;

string_omega_2 = strrep(string_omega_2,'dvarphi_1','result(:,4)') ;
string_omega_2 = strrep(string_omega_2,'varphi_1','result(:,2)') ;
string_omega_2 = strrep(string_omega_2,'dr_1','result(:,3)') ;
string_omega_2 = strrep(string_omega_2,'r_1','result(:,1)') ;
string_omega_2 = strrep(string_omega_2,'^','.^') ;
string_omega_2 = strrep(string_omega_2,'*','.*') ;
string_omega_2 = strrep(string_omega_2,'/','./') ;


x_1=eval(string_x_1);   y_1=eval(string_y_1);   %乙球坐标
x_2=eval(string_x_2);   y_2=eval(string_y_2);   %甲球坐标
omega_1=eval(string_omega_1);   %乙球角速度
omega_2=eval(string_omega_2);   %甲球角速度
theta=eval(string_theta);       %弹簧张角

function dy=ode(t,y)
global s1 s2;
dy=zeros(4,1);
dy(1)=y(3);
dy(2)=y(4);
dy(3)=eval(s1);
dy(4)=eval(s2);
end
%导入微分方程，通过全局变量s1，s2




function fun=dt_diff(fun0)
syms r_1 varphi_1 t;
syms dr_1 dvarphi_1;
syms ddr_1 ddvarphi_1;
fun=diff(fun0,r_1)*dr_1+diff(fun0,dr_1)*ddr_1+...
    diff(fun0,varphi_1)*dvarphi_1+diff(fun0,dvarphi_1)*ddvarphi_1+...
    diff(fun0,t);
end
%求t的全导数函数



