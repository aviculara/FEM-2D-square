function [Up] = assemble_displacements(ENL,NL)
NoN = size(NL,1);
PD = size(NL,2);
DOCs = 0;

for i = 1:NoN
    
    for j = 1:PD
        if (ENL(i,PD+j) == -1)
            DOCs = DOCs + 1;
            Up(DOCs,1) = ENL(i,4*PD+j);
        end
    end
end

end

