function output=testfun2(A11,A12,A21,A22,b11,b12,b21,b22,x1,x2,alpha12,alpha22)

s=1;
d11=s*(A11(1,:)*[x1;x2]+b11(1));
d12=s*alpha12*(A12(1,:)*[x1;x2]+b12(1));
d21=s*(A21(2,:)*[x1;x2]+b21(2));
d22=s*alpha22*(A22(2,:)*[x1;x2]+b22(2));
if d11*d12>0 && abs(d11)>abs(d12)
    f1=d12;
elseif d11*d12>0 && abs(d11)<=abs(d12)
    f1=d11;
else
    f1=0;
end
if d21*d22>0 && abs(d21)>abs(d22)
    f2=d22;
elseif d21*d22>0 && abs(d21)<=abs(d22)
    f2=d21;
else
    f2=0;
end
output=[f1;f2];
end