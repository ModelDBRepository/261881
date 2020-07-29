function g=impl_inf_deconv(f,h)

l=min(length(f),length(h));
A=zeros(l,l); b=zeros(l,1);
A(1,1)=f(1);
for eq=2:l
    for i=0:eq-1
        A(eq,i+1)=f(eq-i);
    end
end
b=h';
g=(A\b)';