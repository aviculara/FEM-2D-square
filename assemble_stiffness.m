function [K] = assemble_stiffness( ENL ,EL , NL , p , m , inclusion )

NoE = size(EL,1);
NPE = size(EL,2);
if (NPE == 4) || (NPE == 9) || (NPE == 8)
    NoE_void = 4*p*m;
else
    NoE_void = 8*p*m;
end

NoN = size(NL,1);
PD = size(NL,2);

K = zeros(NoN*PD, NoN*PD); 
k=[];
for i = 1:NoE
    nl = EL(i,1:NPE); 
    element_number = i;
    k(i,:) = element_stiffness(nl,NL,element_number,NoE_void,inclusion);
end
    row = ENL(nl(1:NPE), (1:PD)+3*PD);
    column = ENL(nl(1:NPE), (1:PD)+3*PD);
    value = k( 1:NoE,((1:NPE)-1)*PD + (1:PD), ((1:NPE)-1)*PD + (1:PD) );
    K(row,column) = K(row,column) + value;
end

