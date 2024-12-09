A11=[-4,1;1,-10];
b11_temp=5;
A12=[-4,-1;-1,-3];
b12_temp=30;

A21=8*[-2,1;1,-2];
b21_temp=0;
A22=[-5,-1;-1,-2];
b22_temp=0;
x11_temp=20;
x12_temp=35;
x21_temp=30;
x22_temp=40;

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


figure;
xa=-10:40;
y11=(-A11(1,1)*xa-b11_temp)/A11(1,2);
y12=(-A12(1,1)*xa-b12_temp)/A12(1,2);
y21=(-A21(2,1)*xa-b21_temp)/A21(2,2);
y22=(-A22(2,1)*xa-b22_temp)/A22(2,2);
plot(xa,y11,'r--',xa,y12,'r',xa,y21,'b--',xa,y22,'b');

xlim([0,50]);
ylim([-120,80]);
hold on
z1=[];
z2=[];
for w11=0:0.05:1
    for w12=0:0.05:(1-w11)
        for w21=0:0.05:(1-w11-w12)
            for w22=0:0.05:(1-w11-w12-w21)
            A=w11*A11+w12*A12+w21*A21+w22*A22;
            b=w11*b11+w12*b12+w21*b21+w22*b22;
            ztemp=-A^(-1)*b;
            z1=[z1,ztemp(1)];
            z2=[z2,ztemp(2)];
            end
        end
    end
end

% for w11=0
%     for w12=0
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
scatter(z1,z2,'y');

%{
%level set for J21
xp21=[x21(1)-15:0.5:x21(1)];
yp21=[x21(2)-20:0.5:x21(2)];
[X21,Y21] = meshgrid(xp21,yp21);
J21=A21(1,1)/2*X21.^2+A21(2,2)/2*Y21.^2+A21(1,2)*X21.*Y21+b21(1)*X21+b21(2)*Y21;
contour(X21,Y21,J21,15,'c:')

%level set for J22
xp22=[x22(1)-5:0.5:x22(1)+5];
yp22=[x22(2)-10:0.5:x22(2)+10];
[X22,Y22] = meshgrid(xp22,yp22);
J22=A22(1,1)/2*X22.^2+A22(2,2)/2*Y22.^2+A22(1,2)*X22.*Y22+b22(1)*X22+b22(2)*Y22;
contour(X22,Y22,J22,15,'b:')

%level set for J11
xp22=[x11(1)-5:0.5:x11(1)+5];
yp22=[x11(2)-10:0.5:x11(2)+10];
[X22,Y22] = meshgrid(xp22,yp22);
J11=A11(1,1)/2*X22.^2+A11(2,2)/2*Y22.^2+A11(1,2)*X22.*Y22+b11(1)*X22+b11(2)*Y22;
contour(X22,Y22,J11,15,'m--')

%level set for J12
xp22=[x12(1)-10:0.5:x12(1)+10];
yp22=[x12(2)-5:0.5:x12(2)+10];
[X22,Y22] = meshgrid(xp22,yp22);
J12=A12(1,1)/2*X22.^2+A12(2,2)/2*Y22.^2+A12(1,2)*X22.*Y22+b12(1)*X22+b12(2)*Y22;
[C,h]=contour(X22,Y22,J12,15,'--');
h.EdgeColor='r';
%}

legend('${\rm BR}^1_1$','${\rm BR}^1_2$','${\rm BR}^2_1$','${\rm BR}^2_2$','Pareto set');
hold off