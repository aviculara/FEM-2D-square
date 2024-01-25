function [Fp] = assemble_forces( ENL , NL) 

NoN = size(NL,1);
PD = size(NL,2);
DOFs = 0;

for i = 1:NoN
    
    for j = 1:PD
    
        if (ENL(i,PD+j) == 1)
            DOFs = DOFs + 1;
            Fp(DOFs,1) = ENL(i,5*PD+j);
        end
    end
end

end

