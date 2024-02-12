clc;clear;
beta=zeros(1,8);L=zeros(1,8);%设定beta角和每次行走后到原点的距离L，L(1)=0;
cosbeta=zeros(1,8);sinbeta=zeros(1,8);%记录sin，cosbeta 的值，避免重复运算。
x_0=zeros(8,8);y_0=zeros(8,8);z_0=zeros(8,8);%船的位置在水底的垂直投影，x0,y0是位置，z0是深度
theta_half=60;%扫描角theta的一半60deg
costhe=cosd(theta_half);sinthe=sind(theta_half);% theta的sin 和cos值
alph=1.5;%海底倾斜角alpha
tanalph=tand(alph);%tan alpha的值

for i=2:8 %计算beta角与L预处理
    beta(i)=beta(i-1)+45;
    L(i)=1852*0.3+L(i-1);
end
for i=1:8 %i是第8个不同的beta角，0,45，90...315
    cosbeta(i)=cosd(beta(i));sinbeta(i)=sind(beta(i));
    phir=beta(i)-90;phil=beta(i)+90; %沿着船行走方向的左边和右边90度处与x轴夹角。
    x_0(i,1)=0;y_0(i,1)=0;z_0(i,1)=120;%船位置初始化
    sinphil=sind(phil);sinphir=sind(phir);
    cosphil=cosd(phil);cosphir=cosd(phir);
    for j=1:8 %地i个beta角中的8个点：0,0.3,0.6...，2.1
        x_0(i,j)=L(j)*cosbeta(i);%船位置x0y0计算
        y_0(i,j)=L(j)*sinbeta(i);
        z_0(i,j)=z_0(i,1)+x_0(i,j)*tanalph;%水深z0
        tl=z_0(i,j)/(costhe-cosphil*sinthe*tanalph);%船行走方向左边三角形斜边长度
        Wl=sqrt(z_0(i,j)^2+tl^2-2*tl*z_0(i,j)*costhe);%余弦定理
        tr=z_0(i,j)/(costhe-cosphir*sinthe*tanalph);%右边三角形斜边长
        Wr=sqrt(z_0(i,j)^2+tr^2-2*tr*z_0(i,j)*costhe);
        W(i,j)=Wl+Wr;
    end
end
