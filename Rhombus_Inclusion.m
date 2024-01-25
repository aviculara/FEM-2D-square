function [ NL,K ] = Rhombus_Inclusion(d1,d2,p,m,R, element_type,NL_old,EL_old)

    t = 0.1*R;
    corner1 = [(d1/2)-t , (d2/2)-t];
    corner2 = [(d1/2)+t , (d2/2)-t];
    corner3 = [(d1/2)-t , (d2/2)+t];
    corner4 = [(d1/2)+t , (d2/2)+t];
    q = [corner1; corner2; corner3; corner4 ];
    
    PD = 2;
    
    NoN_old = size(NL_old,1);
    
    NoN = (2*(p+1)*(m+1) + 2*(p-1)*(m+1)) - 4*p ;
    
    NoE = 4*p*m;
    
    NPE = 4; %D2QU4N
    
    NL = zeros(NoN,PD);
    
    a = (q(2,1)-q(1,1))/p; 
    b = (q(3,2)-q(1,2))/p; 
    
    coor11 = zeros((p+1)*(m+1),PD) ;
    
    for i = 1:p+1
        
        coor11(i,1) = q(1,1) + (i-1)*a;
        coor11(i,2) = q(1,2);
    end
    
    for i = 1:p+1

        coor11(m*(p+1)+i,1) = (R/2)* (i-1)/p + (d1/2)-(R/4); 
        coor11(m*(p+1)+i,2) =(R/2)*(abs((i)-(p+2)/2)/(p))+(d2/2)-R/2 ;
    
     end 
        
      for i=1:m-1
          
          for j=1:p+1
              
              dx = ( coor11(m*(p+1)+j,1 ) - coor11(j,1) )/m ;
              dy = ( coor11(m*(p+1)+j,2 ) - coor11(j,2) )/m ;
              
              coor11(i*(p+1)+j,1) = coor11((i-1)*(p+1)+j,1) + dx ; 
              coor11(i*(p+1)+j,2) = coor11((i-1)*(p+1)+j,2) + dy ;
              
          end
      end

    
    coor22 = zeros((p+1)*(m+1),PD) ; 
    
    for i = 1:p+1
        
        coor22(i,1) = q(3,1) + (i-1)*a;
        coor22(i,2) = q(3,2);
    end
    
    
    for i = 1:p+1

        coor22(m*(p+1)+i,1) = (R/2)* (i-1)/p + (d1/2)-(R/4);
        coor22(m*(p+1)+i,2) = (R/4)*(1-2*(abs((i)-(p+2)/2)/(p)))+(d2/2)+R/4  ; 
    
     end 
        
      for i=1:m-1 
          
          for j=1:p+1
              
              dx = ( coor22(m*(p+1)+j,1 ) - coor22(j,1) )/m ;
              dy = ( coor22(m*(p+1)+j,2 ) - coor22(j,2) )/m ;
              
              coor22(i*(p+1)+j,1) = coor22((i-1)*(p+1)+j,1) + dx ;
              coor22(i*(p+1)+j,2) = coor22((i-1)*(p+1)+j,2) + dy ;
              
          end
      end
    
    coor33 = zeros((p-1)*(m+1),PD) ;
    
    for i = 1:p-1 
        
        coor33(i,1) = q(1,1); 
        coor33(i,2) = q(1,2) + i*b ;
    end
    
    for i = 1:p-1

        coor33(m*(p-1)+i,1) = (R/4)*(2*(abs((i)-(p)/2)/(p)))+(d1/2)-R/2 ; 
        coor33(m*(p-1)+i,2) = (R/2)* (i)/p + (d2/2)-(R/4); 
    
     end 
        
      for i=1:m-1 
          
          for j=1:p-1
              
              dx = ( coor33(m*(p-1)+j,1 ) - coor33(j,1) )/m ; %:No difference, but dy is negative this time [curve-straight line] :%
              dy = ( coor33(m*(p-1)+j,2 ) - coor33(j,2) )/m ;
              
              coor33(i*(p-1)+j,1) = coor33((i-1)*(p-1)+j,1) + dx ; %:Reading info from one line below :%
              coor33(i*(p-1)+j,2) = coor33((i-1)*(p-1)+j,2) + dy ;
              
          end
      end
    
    coor44 = zeros((p-1)*(m+1),PD) ;
    
    for i = 1:p-1 
        
        coor44(i,1) = q(2,1); 
        coor44(i,2) = q(2,2) + i*b ;

    end
    
    
    for i = 1:p-1

        coor44(m*(p-1)+i,1) = (R/4)*(1-2*(abs((i)-(p)/2)/(p)))+(d1/2)+R/4  ; 
        coor44(m*(p-1)+i,2) = (R/2)* (i)/p + (d2/2)-(R/4) ; 
    
     end 
        
      for i=1:m-1 
          
          for j=1:p-1
              
              dx = ( coor44(m*(p-1)+j,1 ) - coor44(j,1) )/m ; 
              dy = ( coor44(m*(p-1)+j,2 ) - coor44(j,2) )/m ;
              
              coor44(i*(p-1)+j,1) = coor44((i-1)*(p-1)+j,1) + dx ; 
              coor44(i*(p-1)+j,2) = coor44((i-1)*(p-1)+j,2) + dy ;
              
          end
      end
    
   for i=1:m+1
       
       NL((i-1)*4*p+1:i*4*p,:) = [  coor11((i-1)*(p+1)+1:i*(p+1),:);
                                    coor44((i-1)*(p-1)+1:i*(p-1),:);
                                    flipud(coor22((i-1)*(p+1)+1:i*(p+1),:));
                                    flipud(coor33((i-1)*(p-1)+1:i*(p-1),:))];
   end
   
   NL((1+(4*p*m)):((4*p)+(4*p*m)),:) = [];
end


