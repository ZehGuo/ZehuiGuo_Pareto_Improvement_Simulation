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

figure;
xa=[-10:40];
y11=(-A11(1,1)*xa-b11_temp)/A11(1,2);
y12=(-A12(1,1)*xa-b12_temp)/A12(1,2);
y21=(-A21(2,1)*xa-b21_temp)/A21(2,2);
y22=(-A22(2,1)*xa-b22_temp)/A22(2,2);
Uy1=(-UA11(1,1)*xa-Ub11(1))/UA11(1,2);
Uy2=(-UA21(2,1)*xa-Ub21(2))/UA21(2,2);
plot(xa,y11,'r--',xa,y12,'r',xa,y21,'b--',xa,y22,'b');
xlim([-5,30]);
ylim([-15,15]);

hold on

%scatter([x11(1),x12(1),x21(1),x22(1)],[x11(2),x12(2),x21(2),x22(2)]);
scatter(ztemp(1),ztemp(2));
fimplicit(@(x,y) ([x,y]*UA*[x;y]/2+Ub'*[x;y])-(x0'*UA*x0/2+Ub'*x0),[-50 50 -50 50],'m--');
fimplicit(@(x,y) ([x,y]*(A11+A21)*[x;y]/2+(b11'+b21')*[x;y])-(x0'*(A11+A21)*x0/2+(b11'+b21')*x0),[-50 50 -50 50],'m');
fimplicit(@(x,y) ([x,y]*(A11+A21)*[x;y]/2+(b11'+b21')*[x;y]-x0'*(A11+A21)*x0/2-(b11'+b21')*x0)-([x,y]*(UA11+UA21)*[x;y]/2+(Ub11'+Ub21')*[x;y]-x0'*UA*x0/2-Ub'*x0),[-50 50 -50 50],'m:');

legend('${\rm BR}^1_1$','${\rm BR}^1_2$','${\rm BR}^2_1$','${\rm BR}^2_2$','$\hat{x}^*$','$\mathcal{D}(\tilde{U},x_0)$','$\mathcal{D}(J_{\rm sum},x_0)$','$\mathcal{D}_{\rm bud}(\sigma)$');
hold off