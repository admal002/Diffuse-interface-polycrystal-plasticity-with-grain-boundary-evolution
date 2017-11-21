function theta=latticeRotation(F11,F12,F21,F22)
theta=zeros(length(F11),1);
for n=1:length(F11)
    F=[F11(n) F12(n);F21(n) F22(n)];
    [R U V] = poldecomp(F);
    theta(n)=atan2(R(2,1),R(1,1));
end