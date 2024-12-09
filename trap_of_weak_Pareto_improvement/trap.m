A11=[-2,1;1,-3];
b11_temp=5;
A12=[-2,-1;-1,-10];
b12_temp=20;
A21=[-4,1;1,-4];
%A21=[-5,1;1,-4];
b21_temp=0;
A22=[-5,-1;-1,-2];
b22_temp=0;
x11_temp=2;
x12_temp=4;
x21_temp=24;
x22_temp=16;


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
%max point for sum J11+J21
xp1=-inv(A11+A21)*(b11+b21);



alpha2=1;
t0=0;
tfinal=4;
%x0=[17;-8.6];
%x0=[18;-8.7];
%x0=[18.5;-8.7];
x0=[18.5;-8.71];

%x0=[20;-10];
[t,xt]=ode45(@(t,xt) testfun2(A11,A12,A21,A22,b11,b12,b21,b22,xt(1),xt(2),1.2*A11(1,1)/A12(1,1),1*A21(2,2)/A22(2,2)),[t0 tfinal],x0);

t1=0.5;
t2=3.35;
[~, idx] = min(abs(t - t2));

figure;
xa=[-10:30];
y11=(-A11(1,1)*xa-b11_temp)/A11(1,2);
y12=(-A12(1,1)*xa-b12_temp)/A12(1,2);
y21=(-A21(2,1)*xa-b21_temp)/A21(2,2);
y22=(-A22(2,1)*xa-b22_temp)/A22(2,2);
plot(xa,y11,'m',xa,y12,'r',xa,y21,'k',xa,y22,'b');
xlim([12,22]);
ylim([-12,-2]);
% xlim([0,30]);
% ylim([-15,15]);


hold on

%level set for J21
% xp22=[x21(1)-5:0.5:x21(1)+5];
% yp22=[x21(2)-10:0.5:x21(2)+10];
xp22=[12:0.1:22];
yp22=[-12:0.1:-2];
[X22,Y22] = meshgrid(xp22,yp22);
J21=A21(1,1)/2*X22.^2+A21(2,2)/2*Y22.^2+A21(1,2)*X22.*Y22+b21(1)*X22+b21(2)*Y22;
contour(X22,Y22,J21,15,'k:')

%level set for J22
% xp22=[x22(1)-5:0.5:x22(1)+5];
% yp22=[x22(2)-10:0.5:x22(2)+10];
xp22=[12:0.1:22];
yp22=[-12:0.1:-2];
[X22,Y22] = meshgrid(xp22,yp22);
J22=A22(1,1)/2*X22.^2+A22(2,2)/2*Y22.^2+A22(1,2)*X22.*Y22+b22(1)*X22+b22(2)*Y22;
contour(X22,Y22,J22,15,'b:')

%scatter([x0(1),xt(t==t1,1),xt(idx,1)],[x0(2),xt(t==t1,2),xt(idx,2)]);


fimplicit(@(x,y) (A21(1,:)*[x;y]+b21(1))*testfun2_1(A11,A12,b11,b12,x,y,1.2*A11(1,1)/A12(1,1))+alpha2*(A21(2,:)*[x;y]+b21(2))*testfun2_2(A21,A22,b21,b22,x,y,A21(2,2)/A22(2,2)),[-30 50 -50 40],'k--');
fimplicit(@(x,y) (A22(1,:)*[x;y]+b22(1))*testfun2_1(A11,A12,b11,b12,x,y,1.2*A11(1,1)/A12(1,1))+alpha2*(A22(2,:)*[x;y]+b22(2))*testfun2_2(A21,A22,b21,b22,x,y,A21(2,2)/A22(2,2)),[-30 50 -50 40],'b--');

plot(xt(:,1),xt(:,2));

%legend('$\hat{\rm BR}^{11}$','$\hat{\rm BR}^{12}$','$\hat{\rm BR}^{21}$','$\hat{\rm BR}^{22}$','$\hat{J}^2_1$','$\hat{J}^2_2$','$(\hat{\Omega}^{21})^{\mathsf{c}}$','$(\hat{\Omega}^{22})^{\mathsf{c}}$','trajectory from $x_0$');
hold off



figure;
xt1=xt(:,1);
xt2=xt(:,2);
j210=x0'*A21*x0/2+b21'*x0;
j220=x0'*A22*x0/2+b22'*x0;
j21=A21(1,1)/2*xt1.^2+A21(2,2)/2*xt2.^2+A21(1,2)*xt1.*xt2+b21(1)*xt1+b21(2)*xt2-x0'*A21*x0/2-b21'*x0;
j22=A22(1,1)/2*xt1.^2+A22(2,2)/2*xt2.^2+A22(1,2)*xt1.*xt2+b22(1)*xt1+b22(2)*xt2-x0'*A22*x0/2-b22'*x0;
plot(t,j21,t,j22);
% hold on
% xline(t1);
% xline(t(idx));
hold off

% figure;
% j110=x0'*A11*x0/2+b11'*x0;
% j120=x0'*A12*x0/2+b12'*x0;
% j11=A11(1,1)/2*xt1.^2+A11(2,2)/2*xt2.^2+A11(1,2)*xt1.*xt2+b11(1)*xt1+b11(2)*xt2-x0'*A11*x0/2-b11'*x0;
% j12=A12(1,1)/2*xt1.^2+A12(2,2)/2*xt2.^2+A12(1,2)*xt1.*xt2+b12(1)*xt1+b12(2)*xt2-x0'*A12*x0/2-b12'*x0;
% plot(t,j11,t,j12);
