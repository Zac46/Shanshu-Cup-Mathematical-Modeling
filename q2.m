clc;clear;
beta=zeros(1,8);L=zeros(1,8);%�趨beta�Ǻ�ÿ�����ߺ�ԭ��ľ���L��L(1)=0;
cosbeta=zeros(1,8);sinbeta=zeros(1,8);%��¼sin��cosbeta ��ֵ�������ظ����㡣
x_0=zeros(8,8);y_0=zeros(8,8);z_0=zeros(8,8);%����λ����ˮ�׵Ĵ�ֱͶӰ��x0,y0��λ�ã�z0�����
theta_half=60;%ɨ���theta��һ��60deg
costhe=cosd(theta_half);sinthe=sind(theta_half);% theta��sin ��cosֵ
alph=1.5;%������б��alpha
tanalph=tand(alph);%tan alpha��ֵ

for i=2:8 %����beta����LԤ����
    beta(i)=beta(i-1)+45;
    L(i)=1852*0.3+L(i-1);
end
for i=1:8 %i�ǵ�8����ͬ��beta�ǣ�0,45��90...315
    cosbeta(i)=cosd(beta(i));sinbeta(i)=sind(beta(i));
    phir=beta(i)-90;phil=beta(i)+90; %���Ŵ����߷������ߺ��ұ�90�ȴ���x��нǡ�
    x_0(i,1)=0;y_0(i,1)=0;z_0(i,1)=120;%��λ�ó�ʼ��
    sinphil=sind(phil);sinphir=sind(phir);
    cosphil=cosd(phil);cosphir=cosd(phir);
    for j=1:8 %��i��beta���е�8���㣺0,0.3,0.6...��2.1
        x_0(i,j)=L(j)*cosbeta(i);%��λ��x0y0����
        y_0(i,j)=L(j)*sinbeta(i);
        z_0(i,j)=z_0(i,1)+x_0(i,j)*tanalph;%ˮ��z0
        tl=z_0(i,j)/(costhe-cosphil*sinthe*tanalph);%�����߷������������б�߳���
        Wl=sqrt(z_0(i,j)^2+tl^2-2*tl*z_0(i,j)*costhe);%���Ҷ���
        tr=z_0(i,j)/(costhe-cosphir*sinthe*tanalph);%�ұ�������б�߳�
        Wr=sqrt(z_0(i,j)^2+tr^2-2*tr*z_0(i,j)*costhe);
        W(i,j)=Wl+Wr;
    end
end
