function f2=testfun2_2(A21,A22,b21,b22,x1,x2,alpha22)

d21=1*(A21(2,:)*[x1;x2]+b21(2));
d22=alpha22*(A22(2,:)*[x1;x2]+b22(2));

if d21*d22>0 && abs(d21)>abs(d22)
    f2=d22;
elseif d21*d22>0 && abs(d21)<=abs(d22)
    f2=d21;
else
    f2=0;
end

end