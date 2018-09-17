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
                    
                  %generate combinations of periodic boundaries for 3d, to check if point
                 %belongs in the interior of any ellipse
                 optionArrX =[i, i+64, i-64];
                 optionArrY =[j, j+64, j-64];
                 optionArrZ =[k, k+64, k-64];

                 ctr = 1;
                 for xi = 1:3
                     for yi=1:3
                         for zi=1:3
                             pts_to_check(ctr,:) = [optionArrX(xi), optionArrY(yi), optionArrZ(zi)];
                             ctr = ctr+ 1;
                         end
                     end
                 end
                    [result, distance] = checkIfEllipseGlobal(pts_to_check,x0,y0,z0,a,b,c,psi1,psi2,phi);
                    
                    out(i,j,k)  = out(i,j,k) + result; 
                    dista(i,j,k) = distance;
                end
            end
        end
        
        %still have to account for periodic boundary condition
        
   end

 
   %% visualize 
   
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

