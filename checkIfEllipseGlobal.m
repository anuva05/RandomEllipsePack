function [out, distance]= checkIfEllipseGlobal(pts_to_check,x0, y0, z0, a,b,c, psi1, psi2, phi)

thresh = 1;


% set diameter and radius
if a > b; diameter = 2*a; else diameter = 2*b; end
if 2*c > diameter; diameter = 2*c; end
if mod(diameter,round(diameter))~= 0; diameter = ceil(diameter); end
radius = diameter/2;
if mod(radius,round(radius))~= 0;
    diameter = diameter+1; radius = diameter/2;
end




% M is the rotation matrix for the ellipse, based on the Euler angles
% psi1, psi2, and phi.  Different rotation matrices and Euler angle
% conventions can be used, if necessary.
M = zeros(3);
M(1,1) = cos(psi1)*cos(psi2)-sin(psi1)*sin(psi2)*cos(phi);
M(2,1) = -cos(psi1)*sin(psi2)-sin(psi1)*cos(psi2)*cos(phi);
M(3,1) = sin(psi1)*sin(phi);
M(1,2) = sin(psi1)*cos(psi2)+cos(psi1)*sin(psi2)*cos(phi);
M(2,2) = -sin(psi1)*sin(psi2)+cos(psi1)*cos(psi2)*cos(phi);
M(3,2) = -cos(psi1)*sin(phi);
M(1,3) = sin(psi2)*sin(phi);
M(2,3) = cos(psi2)*sin(phi);
M(3,3) = cos(phi);




 %% for periodic boundary, do #1



 
%  
%    pts_to_check =[i j k;
%      i+64, j, k; 
%      i, j + 64, k;
%      i, j , k+64;
%      i+64, j+64, k;
%      i+64, j, k+64;
%      i, j+64,k+64;
%      i+64,j+64,k+64;
%      i-64, j, k; 
%      i, j - 64, k ;
%     i, j , k - 64;
%     i-64, j-64, k;
%     i-64, j , k-64;
%     i, j-64,k-64;
%     i-64,j-64,k-64;
%     i+64,j+64,k-64;
%     i+64, j-64,k-64; 
%     i-64,j+64,k-64;
%     i-64,j+64,k+64;
%     i-64,j-64,k+64;
%     i+64,j-64,k+64;
%     i+64; j] ;



%there are more cases... +64 and - 64 together..

    
    
for u = 1:size(pts_to_check,1)
    ptx = pts_to_check(u,1);
    pty = pts_to_check(u,2);
    ptz = pts_to_check(u,3);
    
    
    % distance of given point from centroid of ellipse in question
    a1 = [ptx-x0, pty-y0, ptz-z0];
    a2 = [a^2, b^2, c^2];
    c1 = a1*M;
    c2 = c1.^2./a2;
    distance = sum(c2);
      
   % if distance is less than or equal to 1, then point is in interior or
   % ellipsoid
   out = distance<=  thresh;
   
   if out>=1
       break 
   end
end

