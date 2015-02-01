mode(0)
clear all

A=[1/sqrt(2);1/sqrt(2)]
B=[1,1;1,-1]
C=[sqrt(3)/sqrt(2),0,-sqrt(3)/sqrt(2);0,sqrt(3)/sqrt(2),-sqrt(3)/sqrt(2)]


M=zeros(3,2)

for i=1:2
    for j=1:2
        B(i,j)=A(i)*B(i,j)
    end
end

a=0
k=1

for cpt=1:2
    for j=1:3
        for i=1:2
            M(j,cpt) = M(j,cpt) + B(cpt,i) * C(i,k)
        end
        k=k+1
    end
    k=1
end

C=[0;0]

for cpt=1:2
    for j=1:3
        C(cpt)= M(j,cpt)^2 + C(cpt)
    end
    C(cpt) = sqrt(C(cpt))
end




