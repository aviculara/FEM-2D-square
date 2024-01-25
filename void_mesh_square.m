function [ NL , EL ] = void_mesh_square(d1,d2,p,m,R, element_type,inclusion)

    q = [0 0; d1 0; 0 d2; d1 d2 ]; 
    
    PD = 2;
    
    NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1) ;
    
    NoE = 4*p*m;
    
    NPE = 4; 
    
    NL = zeros(NoN,PD);
    
    a = (q(2,1)-q(1,1))/p; 
    b = (q(3,2)-q(1,2))/p; 
    
    
    
    coor11 = zeros((p+1)*(m+1),PD) ; 
    
    for i = 1:p+1
        
        coor11(i,1) = q(1,1) + (i-1)*a;
        coor11(i,2) = q(1,2);
    end
    
    
    for i = 1:p+1
        
        coor11(m*(p+1)+i,1) = (d1-R)/2 + (i-1)*( R/p ) ;
        coor11(m*(p+1)+i,2) = (d2-R)/2;
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
        
        coor22(m*(p+1)+i,1) = (d1-R)/2 + (i-1)*( R/p );
        coor22(m*(p+1)+i,2) = d2-(d2-R)/2; 
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
        
        coor33(m*(p-1)+i,1) = (d1-R)/2 ; 
        coor33(m*(p-1)+i,2) = (d2-R)/2 + (i)*(R/p); 

     end 
        
      for i=1:m-1 
          
          for j=1:p-1
              
              dx = ( coor33(m*(p-1)+j,1 ) - coor33(j,1) )/m ; 
              dy = ( coor33(m*(p-1)+j,2 ) - coor33(j,2) )/m ;
              
              coor33(i*(p-1)+j,1) = coor33((i-1)*(p-1)+j,1) + dx ; 
              coor33(i*(p-1)+j,2) = coor33((i-1)*(p-1)+j,2) + dy ;
              
          end
          
      end
    
   
    
    coor44 = zeros((p-1)*(m+1),PD) ;
    
    for i = 1:p-1 
        
        coor44(i,1) = q(2,1); 
        coor44(i,2) = q(2,2) + i*b ;
    end
    
    
    for i = 1:p-1
        
        coor44(m*(p-1)+i,1) = d1-(d1-R)/2 ; 
        coor44(m*(p-1)+i,2) = (d2-R)/2 + (i)*(R/p); 
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
   
   
   
   
   EL = zeros(NoE,NPE);
   
   for i=1:m
       
       for j=1:4*p  
           
           if (j==1)
               
                EL((i-1)*(4*p)+j,1) = (i-1)*(4*p) + 1;
                EL((i-1)*(4*p)+j,2) =  EL((i-1)*(4*p)+j,1) + 1;
                EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p;
                EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1; 
               
           elseif (j==4*p) 
                
                EL((i-1)*(4*p)+j,1) = (i)*(4*p); 
                EL((i-1)*(4*p)+j,2) =  (i-1)*(4*p)+1;
                EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,1) + 1;
                EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p;
               
               
           else
                EL((i-1)*(4*p)+j,1) =  EL((i-1)*(4*p)+j-1,2); 
                EL((i-1)*(4*p)+j,4) =  EL((i-1)*(4*p)+j-1,3);
                EL((i-1)*(4*p)+j,3) =  EL((i-1)*(4*p)+j,4) + 1;
                EL((i-1)*(4*p)+j,2) =  EL((i-1)*(4*p)+j,1) + 1;
                
           end

       end
   end


    
    
    
 
if inclusion == 1
    old_NL = NL;
    old_NoN = size(old_NL,1);
    old_NoN2 = old_NoN -(4*p);
    old_EL = EL;
    [NL] = square_inclusion(d1,d2,p,p,R, element_type,old_NL,old_EL);
    EL_inc = zeros(p^2,NPE);

    for i = 1:p 

            for j = 1:p

                if (j==1) && (i ==1) 

                    EL_inc( (i-1)*p+j , 1 ) = (i-1) * (p+1)+j; 
                    EL_inc( (i-1)*p+j , 2) = EL_inc( (i-1)*p+j , 1 ) + 1 ;               
                    EL_inc( (i-1)*p+j , 3) = EL_inc( (i-1)*p+j , 1 ) + (4*p) ;
                    EL_inc( (i-1)*p+j , 4) = EL_inc( (i-1)*p+j , 3) -1 ;

                elseif (j==1) && (i~=1) && (i<p)
                    EL_inc( (i-1)*p+j , 2 ) = EL_inc( (i-2)*p+j , 3 ); 
                    EL_inc( (i-1)*p+j , 1) =  EL_inc( (i-2)*p+j , 4 );               
                    EL_inc( (i-1)*p+j , 3) = EL_inc( (i-1)*p+j , 2 ) + (p-1) ;
                    EL_inc( (i-1)*p+j , 4) = EL_inc( (i-1)*p+j , 1) -1 ;
                 elseif (j==1) &&(i == p)
                    EL_inc( (i-1)*p+j , 2 ) = EL_inc( (i-2)*p+j , 3 ); 
                    EL_inc( (i-1)*p+j , 1) =  EL_inc( (i-2)*p+j , 4 );               
                    EL_inc( (i-1)*p+j , 4) = EL_inc( (i-1)*p+j , 1) -1;
                    EL_inc( (i-1)*p+j , 3) = EL_inc( (i-1)*p+j , 4) -1;

                elseif (j==p)&& (i==1)
                    EL_inc( (i-1)*p+j , 2 ) =  (p+1); 
                    EL_inc( (i-1)*p+j , 1) = EL_inc( (i-1)*p+j , 2 ) - 1 ;               
                    EL_inc( (i-1)*p+j , 3) = EL_inc( (i-1)*p+j , 2 ) + 1 ;
                    EL_inc( (i-1)*p+j , 4) = EL_inc( (i-1)*p+j , 1)+ ((4*p)-1) ;
                elseif (j==p) && (i~=1) && (i<p)
                    EL_inc( (i-1)*p+j , 2 ) = EL_inc( (i-2)*p+j , 3 ); 
                    EL_inc( (i-1)*p+j , 1) =  EL_inc( (i-2)*p+j , 4 );               
                    EL_inc( (i-1)*p+j , 4) = EL_inc( (i-1)*p+j , 1 ) + (p-1) ;
                    EL_inc( (i-1)*p+j , 3) = EL_inc( (i-1)*p+j , 2) +1 ;
                elseif (j==p) && (i == p)
                    EL_inc( (i-1)*p+j , 2 ) = EL_inc( (i-2)*p+j , 3 ); 
                    EL_inc( (i-1)*p+j , 1) =  EL_inc( (i-2)*p+j , 4 );               
                    EL_inc( (i-1)*p+j , 3) = EL_inc( (i-1)*p+j , 2) +1;
                    EL_inc( (i-1)*p+j , 4) = EL_inc( (i-1)*p+j , 3) +1;

                elseif (j>1) && (j<p) 
                    
                    if (i~=p)
                    EL_inc( (i-1)*p+j , 1 ) = EL_inc( (i-1)*p+(j-1) , 2 ); 
                    EL_inc( (i-1)*p+j , 4 ) = EL_inc( (i-1)*p+(j-1) , 3 );
                    EL_inc( (i-1)*p+j , 2 ) = EL_inc( (i-1)*p+j , 1 ) + 1; 
                    EL_inc( (i-1)*p+j , 3 ) = EL_inc( (i-1)*p+j , 4 ) + 1;
                    elseif i == p
                    EL_inc( (i-1)*p+j , 1 ) = EL_inc( (i-1)*p+(j-1) , 2 ); 
                    EL_inc( (i-1)*p+j , 4 ) = EL_inc( (i-1)*p+(j-1) , 3 );
                    EL_inc( (i-1)*p+j , 2 ) = EL_inc( (i-1)*p+j , 1 ) + 1; 
                    EL_inc( (i-1)*p+j , 3 ) = EL_inc( (i-1)*p+j , 4 ) - 1;


                    
                    end
                end
            end

    end
    disp (EL_inc)
    EL_inc = EL_inc + (4*p*m);
    disp(4*p*m)
    disp(EL_inc)
    EL = vertcat(EL,EL_inc);
end




if isequal(element_type,'D2QU8N')
    gate = 1;
        for i=1:size(EL,1)
            nodenumcount = size(NL,1);
            nnum1 = EL(i,1);
            nnum2 = EL(i,2);
            nnum3 = EL(i,3);
            nnum4 = EL(i,4); 

            ncoord1 = NL(nnum1,:);
            ncoord2 = NL(nnum2,:);
            ncoord3 = NL(nnum3,:);
            ncoord4 = NL(nnum4,:); 

            addnode1 =  [(ncoord1(1,1)+ncoord2(1,1))/2 , (ncoord1(1,2)+ncoord2(1,2))/2];
            addnode2 =  [(ncoord2(1,1)+ncoord3(1,1))/2 , (ncoord2(1,2)+ncoord3(1,2))/2];
            addnode3 =  [(ncoord3(1,1)+ncoord4(1,1))/2 , (ncoord3(1,2)+ncoord4(1,2))/2];
            addnode4 =  [(ncoord4(1,1)+ncoord1(1,1))/2 , (ncoord4(1,2)+ncoord1(1,2))/2];

   

            if ismember(addnode1,NL,'rows') == 0
                NL = [NL;addnode1];
                nodenumcount = nodenumcount + 1;
                nodeinfo1 = [nodenumcount,addnode1];
            else
                nodeinfo1 = [0,0,0];
            end
            if ismember(addnode2,NL,'rows') == 0
                NL = [NL;addnode2];
                nodenumcount = nodenumcount + 1;
                nodeinfo2 = [nodenumcount,addnode2];
            else
                nodeinfo2 = [0,0,0];
            end
            if ismember(addnode3,NL,'rows') == 0
                NL = [NL;addnode3];
                nodenumcount = nodenumcount + 1;
                nodeinfo3 = [nodenumcount,addnode3];
            else
                nodeinfo3 = [0,0,0];
            end
            if ismember(addnode4,NL,'rows') == 0
                NL = [NL;addnode4];
                nodenumcount = nodenumcount + 1;
                nodeinfo4 = [nodenumcount,addnode4];
            else
                nodeinfo4 = [0,0,0];
            end
   
             if gate == 1

                zeroarray = zeros(size(EL,1),4);
                EL = [EL zeroarray];
                gate = 0;
             end

            if not(nodeinfo1 == [0,0,0])
                EL(i,5) = nodeinfo1(1,1);
            else
                RowIdx = find(ismember(NL, addnode1,'rows'));
                EL(i,5) = RowIdx;  
            end

            if not(nodeinfo2 == [0,0,0])
                EL(i,6) = nodeinfo2(1,1);
            else
                RowIdx = find(ismember(NL, addnode2,'rows'));
                EL(i,6) = RowIdx;  
            end

            if not(nodeinfo3 == [0,0,0])
                EL(i,7) = nodeinfo3(1,1);
            else
                RowIdx = find(ismember(NL, addnode3,'rows'));
                EL(i,7) = RowIdx;  
            end

            if not(nodeinfo4 == [0,0,0])
                EL(i,8) = nodeinfo4(1,1);
            else
                RowIdx = find(ismember(NL, addnode4,'rows'));
                EL(i,8) = RowIdx;  
            end

   


        end
end
if isequal(element_type,'D2QU9N')
    gate = 1;
        for i=1:size(EL,1)
            nodenumcount = size(NL,1);
            nnum1 = EL(i,1);
            nnum2 = EL(i,2);
            nnum3 = EL(i,3);
            nnum4 = EL(i,4);

            ncoord1 = NL(nnum1,:);
            ncoord2 = NL(nnum2,:);
            ncoord3 = NL(nnum3,:);
            ncoord4 = NL(nnum4,:); 

            addnode1 =  [(ncoord1(1,1)+ncoord2(1,1))/2 , (ncoord1(1,2)+ncoord2(1,2))/2];
            addnode2 =  [(ncoord2(1,1)+ncoord3(1,1))/2 , (ncoord2(1,2)+ncoord3(1,2))/2];
            addnode3 =  [(ncoord3(1,1)+ncoord4(1,1))/2 , (ncoord3(1,2)+ncoord4(1,2))/2];
            addnode4 =  [(ncoord4(1,1)+ncoord1(1,1))/2 , (ncoord4(1,2)+ncoord1(1,2))/2];

            addnode5 = [(ncoord1(1,1)+ncoord2(1,1)+ncoord3(1,1)+ncoord4(1,1))/4, ((ncoord1(1,2)+ncoord2(1,2)+ncoord3(1,2)+ncoord4(1,2))/4)];

            if ismember(addnode1,NL,'rows') == 0
                NL = [NL;addnode1];
                nodenumcount = nodenumcount + 1;
                nodeinfo1 = [nodenumcount,addnode1];
            else
                nodeinfo1 = [0,0,0];
            end
            if ismember(addnode2,NL,'rows') == 0
                NL = [NL;addnode2];
                nodenumcount = nodenumcount + 1;
                nodeinfo2 = [nodenumcount,addnode2];
            else
                nodeinfo2 = [0,0,0];
            end
            if ismember(addnode3,NL,'rows') == 0
                NL = [NL;addnode3];
                nodenumcount = nodenumcount + 1;
                nodeinfo3 = [nodenumcount,addnode3];
            else
                nodeinfo3 = [0,0,0];
            end
            if ismember(addnode4,NL,'rows') == 0
                NL = [NL;addnode4];
                nodenumcount = nodenumcount + 1;
                nodeinfo4 = [nodenumcount,addnode4];
            else
                nodeinfo4 = [0,0,0];
            end
             if ismember(addnode5,NL,'rows') == 0
                 NL = [NL;addnode5];
                 nodenumcount = nodenumcount + 1;
                 nodeinfo5 = [nodenumcount,addnode5];
             else
                 nodeinfo5 = [0,0,0];
             end
             zeroarray = zeros(size(EL,1),5);
             if gate == 1

                zeroarray = zeros(size(EL,1),5);
                EL = [EL zeroarray];
                gate = 0;
             end

            if not(nodeinfo1 == [0,0,0])
                EL(i,5) = nodeinfo1(1,1);
            else
                RowIdx = find(ismember(NL, addnode1,'rows'));
                EL(i,5) = RowIdx;  
            end

            if not(nodeinfo2 == [0,0,0])
                EL(i,6) = nodeinfo2(1,1);
            else
                RowIdx = find(ismember(NL, addnode2,'rows'));
                EL(i,6) = RowIdx;  
            end

            if not(nodeinfo3 == [0,0,0])
                EL(i,7) = nodeinfo3(1,1);
            else
                RowIdx = find(ismember(NL, addnode3,'rows'));
                EL(i,7) = RowIdx;  
            end

            if not(nodeinfo4 == [0,0,0])
                EL(i,8) = nodeinfo4(1,1);
            else
                RowIdx = find(ismember(NL, addnode4,'rows'));
                EL(i,8) = RowIdx;  
            end

             if not(nodeinfo5 == [0,0,0])
                 EL(i,9) = nodeinfo5(1,1);
             else
                 RowIdx = find(ismember(NL, addnode2,'rows'));
                 EL(i,9) = RowIdx;  
             end


        end
end
NoE = size(EL,1);

    if isequal(element_type,'D2TR3N')
        
        NPE_new = 3; 
        NoE_new = 2*NoE;
        EL_new = zeros(NoE_new , NPE_new);
        
        for i = 1:NoE
            
           EL_new( 2*(i-1)+1 , 1) = EL(i,1);  
           EL_new( 2*(i-1)+1 , 2) = EL(i,2);
           EL_new( 2*(i-1)+1 , 3) = EL(i,3);
           
           EL_new( 2*(i-1)+2 , 1) = EL(i,1);
           EL_new( 2*(i-1)+2 , 2) = EL(i,3);
           EL_new( 2*(i-1)+2 , 3) = EL(i,4);
           
 
        end
    
        EL = EL_new;

    end
    if isequal(element_type,'D2TR6N')
        
        
        gate = 1;
        for i=1:size(EL,1)
            nodenumcount = size(NL,1);
            nnum1 = EL(i,1);
            nnum2 = EL(i,2);
            nnum3 = EL(i,3);
            nnum4 = EL(i,4); 

            ncoord1 = NL(nnum1,:);
            ncoord2 = NL(nnum2,:);
            ncoord3 = NL(nnum3,:);
            ncoord4 = NL(nnum4,:); 
            addnode1 =  [(ncoord1(1,1)+ncoord2(1,1))/2 , (ncoord1(1,2)+ncoord2(1,2))/2];
            addnode2 =  [(ncoord2(1,1)+ncoord3(1,1))/2 , (ncoord2(1,2)+ncoord3(1,2))/2];
            addnode3 =  [(ncoord3(1,1)+ncoord4(1,1))/2 , (ncoord3(1,2)+ncoord4(1,2))/2];
            addnode4 =  [(ncoord4(1,1)+ncoord1(1,1))/2 , (ncoord4(1,2)+ncoord1(1,2))/2];

   

            if ismember(addnode1,NL,'rows') == 0
                NL = [NL;addnode1];
                nodenumcount = nodenumcount + 1;
                nodeinfo1 = [nodenumcount,addnode1];
            else
                nodeinfo1 = [0,0,0];
            end
            if ismember(addnode2,NL,'rows') == 0
                NL = [NL;addnode2];
                nodenumcount = nodenumcount + 1;
                nodeinfo2 = [nodenumcount,addnode2];
            else
                nodeinfo2 = [0,0,0];
            end
            if ismember(addnode3,NL,'rows') == 0
                NL = [NL;addnode3];
                nodenumcount = nodenumcount + 1;
                nodeinfo3 = [nodenumcount,addnode3];
            else
                nodeinfo3 = [0,0,0];
            end
            if ismember(addnode4,NL,'rows') == 0
                NL = [NL;addnode4];
                nodenumcount = nodenumcount + 1;
                nodeinfo4 = [nodenumcount,addnode4];
            else
                nodeinfo4 = [0,0,0];
            end
   
             if gate == 1

                zeroarray = zeros(size(EL,1),4);
                EL = [EL zeroarray];
                gate = 0;
             end

            if not(nodeinfo1 == [0,0,0])
                EL(i,5) = nodeinfo1(1,1);
            else
                RowIdx = find(ismember(NL, addnode1,'rows'));
                EL(i,5) = RowIdx;  
            end

            if not(nodeinfo2 == [0,0,0])
                EL(i,6) = nodeinfo2(1,1);
            else
                RowIdx = find(ismember(NL, addnode2,'rows'));
                EL(i,6) = RowIdx;  
            end

            if not(nodeinfo3 == [0,0,0])
                EL(i,7) = nodeinfo3(1,1);
            else
                RowIdx = find(ismember(NL, addnode3,'rows'));
                EL(i,7) = RowIdx;  
            end

            if not(nodeinfo4 == [0,0,0])
                EL(i,8) = nodeinfo4(1,1);
            else
                RowIdx = find(ismember(NL, addnode4,'rows'));
                EL(i,8) = RowIdx;  
            end

    


        end
        
      
        
        
        
        
        
        NPE_new = 6; 
        NoE_new = 2*NoE;
        EL_new = zeros(NoE_new , NPE_new);
        
        for i = 1:NoE
           avgnum1 = EL(i,1);
           avgnum2 = EL(i,3);
           
           avgnode1 = NL(avgnum1,:);
           avgnode2 = NL(avgnum2,:);
           
           avgcoord = [(avgnode1(1,1)+avgnode2(1,1))/2 (avgnode1(1,2)+avgnode2(1,2))/2];
           NL = [NL ; avgcoord];
           nodenumber = size(NL,1);
           
           
            
           EL_new( 2*(i-1)+1 , 1) = EL(i,1);  
           EL_new( 2*(i-1)+1 , 2) = EL(i,2);
           EL_new( 2*(i-1)+1 , 3) = EL(i,3);
           EL_new( 2*(i-1)+1 , 4) = nodenumber;  
           EL_new( 2*(i-1)+1 , 5) = EL(i,5);
           EL_new( 2*(i-1)+1 , 6) = EL(i,6);
           
           EL_new( 2*(i-1)+2 , 1) = EL(i,1);
           EL_new( 2*(i-1)+2 , 2) = EL(i,3);
           EL_new( 2*(i-1)+2 , 3) = EL(i,4);
           EL_new( 2*(i-1)+2 , 4) = EL(i,8);
           EL_new( 2*(i-1)+2 , 5) = EL(i,7);
           EL_new( 2*(i-1)+2 , 6) = nodenumber;
           
           nodenumber = nodenumber + 1;
        end
        EL = EL_new;
    end



    


end

