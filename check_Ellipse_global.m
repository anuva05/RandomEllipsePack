   close all

   out = zeros(64,64,64);
   for id = 1:7 % size of bounding box is determined by max(a,b,c) and is 2*max(a,b,c)+ 1
        ellipseParams = ellipse(id, :);
        I_ell = collectIellipse(:,:,:,id);
        x0 = ellipseParams(1);
        y0 = ellipseParams(2);
        z0 = ellipseParams(3);
        
        a =  ellipseParams(4);
        b =  ellipseParams(5);
        c=   ellipseParams(6);
        
        psi1 =  ellipseParams(7);
        psi2 =  ellipseParams(8);
        phi =  ellipseParams(9);
     
        for i=1:64
            for j=1:64
                for k=1:64
                    [result, distance] = checkIfEllipseGlobal(i,j,k,x0,y0,z0,a,b,c,psi1,psi2,phi);
                    
                    out(i,j,k)  = out(i,j,k) + result; 
                    dista(i,j,k) = distance;
                end
            end
        end
        
        %still have to account for periodic boundary condition
        
   end
   
   
   for i = 1:64
        image_view(out(:,:,i));
pause(0.15)
   end
    
   
diff = out -I; 

  for i = 1:64
        image_view(diff(:,:,i));
pause(0.15)
  end
    
   
  
    for i = 1:64
        figure(1)
        J = I(:,:,i);
        image_view(J);
title(['Slice ' num2str(i)])
pause

 figure(2)
        J = out(:,:,i);
        image_view(J);
title(['Slice ' num2str(i)])
pause
    end

