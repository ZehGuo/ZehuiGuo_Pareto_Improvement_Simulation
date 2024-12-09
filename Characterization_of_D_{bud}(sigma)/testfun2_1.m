function f1=testfun2_1(A11,A12,b11,b12,x1,x2,alpha12)

d11=1*(A11(1,:)*[x1;x2]+b11(1));
d12=alpha12*(A12(1,:)*[x1;x2]+b12(1));

if d11*d12>0 && abs(d11)>abs(d12)
    f1=d12;
elseif d11*d12>0 && abs(d11)<=abs(d12)
    f1=d11;
else
    f1=0;
end
end