A11=[-4,1;1,-10];
b11_temp=5;
A12=[-4,-1;-1,-3];
b12_temp=30;

A21=8*[-2,1;1,-2];
b21_temp=0;
A22=[-5,-1;-1,-2];
b22_temp=0;
x11_temp=-10;
x12_temp=10;
x21_temp=30;
x22_temp=13;

[temp1,temp2]=testfun1([A11(1,:),b11_temp],A11(2,:),x11_temp);
x11=[x11_temp,temp1];
b11=[b11_temp;temp2];
[temp1,temp2]=testfun1([A12(1,:),b12_temp],A12(2,:),x12_temp);
x12=[x12_temp,temp1];
b12=[b12_temp;temp2];
[temp1,temp2]=testfun1([A21(2,:),b21_temp],A21(1,:),x21_temp);
x21=[x21_temp,temp1];
b21=[temp2;b21_temp];
[temp1,temp2]=testfun1([A22(2,:),b22_temp],A22(1,:),x22_temp);
x22=[x22_temp,temp1];
b22=[temp2;b22_temp];
xp1=-inv(A11+A21)*(b11+b21);


w11=0.7;
w12=4;
w21=1;
w22=2;
A=w11*A11+w12*A12+w21*A21+w22*A22;
b=w11*b11+w12*b12+w21*b21+w22*b22;
ztemp=-A^(-1)*b;

%x0=[20;-15];
x0=[16;-15];
sigma=0.05;
UA=sigma*A;
Ub=sigma*b;
UA11=0.3*UA;
UA21=0.7*UA;
Ub11=0.3*Ub;
Ub21=0.7*Ub;

t0=0;
tfinal=3;
[t,xt]=ode45(@(t,xt) testfun2(A11,A12,A21,A22,b11,b12,b21,b22,xt(1),xt(2),A11(1,1)/A12(1,1),1*A21(2,2)/A22(2,2)),[t0 tfinal],x0);
[ut,uxt]=ode45(@(ut,uxt) testfun2(UA11,A12,UA21,A22,Ub11,b12,Ub21,b22,uxt(1),uxt(2),UA11(1,1)/A12(1,1),1*UA21(2,2)/A22(2,2)),[t0 tfinal],x0);

figure;
xa=[-10:40];
y11=(-A11(1,1)*xa-b11_temp)/A11(1,2);
y12=(-A12(1,1)*xa-b12_temp)/A12(1,2);
y21=(-A21(2,1)*xa-b21_temp)/A21(2,2);
y22=(-A22(2,1)*xa-b22_temp)/A22(2,2);
Uy1=(-UA11(1,1)*xa-Ub11(1))/UA11(1,2);
Uy2=(-UA21(2,1)*xa-Ub21(2))/UA21(2,2);
plot(xa,y11,'r--',xa,y12,'r',xa,Uy1,'r:',xa,y21,'b--',xa,y22,'b',xa,Uy2,'b:');
xlim([-5,30]);
ylim([-15,15]);

hold on

%scatter([x11(1),x12(1),x21(1),x22(1)],[x11(2),x12(2),x21(2),x22(2)]);
%scatter(ztemp(1),ztemp(2));
%fimplicit(@(x,y) ([x,y]*UA*[x;y]/2+Ub'*[x;y])-(x0'*UA*x0/2+Ub'*x0),[-50 50 -50 50],'r--');
%fimplicit(@(x,y) ([x,y]*(A11+A21)*[x;y]/2+(b11'+b21')*[x;y])-(x0'*(A11+A21)*x0/2+(b11'+b21')*x0),[-50 50 -50 50],'r');
fimplicit(@(x,y) ([x,y]*(A11+A21)*[x;y]/2+(b11'+b21')*[x;y]-x0'*(A11+A21)*x0/2-(b11'+b21')*x0)-([x,y]*(UA11+UA21)*[x;y]/2+(Ub11'+Ub21')*[x;y]-x0'*UA*x0/2-Ub'*x0),[-50 50 -50 50],'m:');

plot(xt(:,1),xt(:,2),'k--');
plot(uxt(:,1),uxt(:,2),'k:');

legend('${\rm BR}^1_1$','${\rm BR}^1_2$','$\tilde{\rm BR}^1_1$','${\rm BR}^2_1$','${\rm BR}^2_2$','$\tilde{\rm BR}^2_1$','$\mathcal{D}_{\rm bud}(\sigma)$','original trajectory','incentivized trajectory');

% z1=[];
% z2=[];
% for w11=0:0.25:1
%     for w12=0:0.25:(1-w11)
%         for w21=0:0.25:(1-w11-w12)
%             for w22=0:0.25:(1-w11-w12-w21)
%             A=w11*A11+w12*A12+w21*A21+w22*A22;
%             b=w11*b11+w12*b12+w21*b21+w22*b22;
%             ztemp=-A^(-1)*b;
%             z1=[z1,ztemp(1)];
%             z2=[z2,ztemp(2)];
%             end
%         end
%     end
% end
% scatter(z1,z2);
hold off


figure;
xt1=xt(:,1);
xt2=xt(:,2);
j11=A11(1,1)/2*xt1.^2+A11(2,2)/2*xt2.^2+A11(1,2)*xt1.*xt2+b11(1)*xt1+b11(2)*xt2-x0'*A11*x0/2-b11'*x0;
j12=A12(1,1)/2*xt1.^2+A12(2,2)/2*xt2.^2+A12(1,2)*xt1.*xt2+b12(1)*xt1+b12(2)*xt2-x0'*A12*x0/2-b12'*x0;
j21=A21(1,1)/2*xt1.^2+A21(2,2)/2*xt2.^2+A21(1,2)*xt1.*xt2+b21(1)*xt1+b21(2)*xt2-x0'*A21*x0/2-b21'*x0;
j22=A22(1,1)/2*xt1.^2+A22(2,2)/2*xt2.^2+A22(1,2)*xt1.*xt2+b22(1)*xt1+b22(2)*xt2-x0'*A22*x0/2-b22'*x0;
plot(t,j11/10,t,j12,t,j21/10,t,j22);
legend('${J}^{1}_1(x)$','${J}^{1}_2(x)$','${J}^{2}_1(x)$','${J}^{2}_2(x)$');
%legend('$\hat{J}^{2}_1(x)$','$\hat{J}^{2}_2(x)$');


figure;
uxt1=uxt(:,1);
uxt2=uxt(:,2);
uj11=UA11(1,1)/2*uxt1.^2+UA11(2,2)/2*uxt2.^2+UA11(1,2)*uxt1.*uxt2+Ub11(1)*uxt1+Ub11(2)*uxt2-x0'*UA11*x0/2-Ub11'*x0;
uj12=A12(1,1)/2*uxt1.^2+A12(2,2)/2*uxt2.^2+A12(1,2)*uxt1.*uxt2+b12(1)*uxt1+b12(2)*uxt2-x0'*A12*x0/2-b12'*x0;
uj21=UA21(1,1)/2*uxt1.^2+UA21(2,2)/2*uxt2.^2+UA21(1,2)*uxt1.*uxt2+Ub21(1)*uxt1+Ub21(2)*uxt2-x0'*UA21*x0/2-Ub21'*x0;
uj22=A22(1,1)/2*uxt1.^2+A22(2,2)/2*uxt2.^2+A22(1,2)*uxt1.*uxt2+b22(1)*uxt1+b22(2)*uxt2-x0'*A22*x0/2-b22'*x0;
p=UA(1,1)/2*uxt1.^2+UA(2,2)/2*uxt2.^2+UA(1,2)*uxt1.*uxt2+Ub(1)*uxt1+Ub(2)*uxt2-x0'*UA*x0/2-Ub'*x0-((A11(1,1)+A21(1,1))/2*uxt1.^2+(A11(2,2)+A21(2,2))/2*uxt2.^2+(A11(1,2)+A21(1,2))*uxt1.*uxt2+(b11(1)+b21(1))*uxt1+(b11(2)+b21(2))*uxt2-x0'*(A11+A21)*x0/2-(b11+b21)'*x0);
plot(ut,uj11/10,ut,uj12,ut,uj21/10,ut,uj22,ut,p/10);
legend('$\hat{J}^{1}_1(x)$','${J}^{1}_2(x)$','$\hat{J}^{2}_1(x)$','${J}^{2}_2(x)$','p');