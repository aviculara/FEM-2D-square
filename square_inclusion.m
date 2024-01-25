function [NL ] = square_inclusion(d1,d2,p,m,R, element_type,old_NL,old_EL)
PD = 2 ;
    q = [ ((d1/2)-(R/2)) ((d2/2)-(R/2));((d1/2)+(R/2)) ((d2/2)-(R/2)); ((d1/2)-(R/2)) ((d2/2)+(R/2)); ((d1/2)+(R/2)) ((d2/2)+(R/2))]; %:: Corners of the Domain::%

%     NoN = (p+1) * (m+1);
    NoN = (p-1)^2;
    old_NoN = size(old_NL,1);
    NoE = (p+1)*m;
    
    NPE = 4; %D2QU4N

    NL = zeros(NoN, PD);
    
    a = (q(2,1)-q(1,1))/p ; % a (x increment), b (y increment) 
    b= (q(3,2)-q(1,2))/m ; 
    
    n =1; %counter
    
    for i = 2:m+1-1 % # of nodes
         
        for j = 2:p %First finish the horizontal line
            
            NL(n,1) = q(1,1) + (j-1)*a ;
            NL(n,2) = q(1,2) + (i-1)*b;
            
            n = n+1;
            
        end
    end
    EL = zeros(NoE,NPE);
    
    for i = 1:m % # of elements
        for j = 1:p
            if (j==1) 
                EL( (i-1)*p+j , 1 ) = (i-1) * (p+1)+j; %p elements, p+1 nodes
                EL( (i-1)*p+j , 2) = EL( (i-1)*p+j , 1 ) + 1 ;               
                EL( (i-1)*p+j , 4) = EL( (i-1)*p+j , 1 ) + (p+1) ;
                EL( (i-1)*p+j , 3) = EL( (i-1)*p+j , 4) +1 ;
                
            else
                
                EL( (i-1)*p+j , 1 ) =EL( (i-1)*p+j-1 , 2 ); 
                EL( (i-1)*p+j , 4) = EL( (i-1)*p+j-1 , 3 );
                
                EL( (i-1)*p+j , 2) = EL( (i-1)*p+j , 1 ) + 1; 
                EL( (i-1)*p+j , 3) = EL( (i-1)*p+j , 4 ) + 1;
            
            end
        end
        
        
    end
%     if isequal(element_type,'D2TR3N')
%         
%         NPE_new = 3; %triangular
%         NoE_new = 2*NoE;
%         EL_new = zeros(NoE_new , NPE_new);
%         
%         for i = 1:NoE
%            EL_new( 2*(i-1)+1 , 1) = EL(i,1);  
%            EL_new( 2*(i-1)+1 , 2) = EL(i,2);
%            EL_new( 2*(i-1)+1 , 3) = EL(i,3);

%            EL_new( 2*(i-1)+2 , 1) = EL(i,1);
%            EL_new( 2*(i-1)+2 , 2) = EL(i,3);
%            EL_new( 2*(i-1)+2 , 3) = EL(i,4);
%         end
%     
%         EL = EL_new;
%     end

EL = EL + old_NoN;
NL = vertcat(old_NL,NL);
EL = vertcat(old_EL,EL);

end

