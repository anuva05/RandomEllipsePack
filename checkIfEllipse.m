function [out, distance]= checkIfEllipse(x0, y0, z0, a,b,c, psi1, psi2, phi)

thresh = 1;


% set diameter and radius
if a > b; diameter = 2*a; else diameter = 2*b; end
if 2*c > diameter; diameter = 2*c; end
if mod(diameter,round(diameter))~= 0; diameter = ceil(diameter); end
radius = diameter/2;
if mod(radius,round(radius))~= 0;
    diameter = diameter+1; radius = diameter/2;
end

%pre-set dist values to 2
%dist(1:diameter+1,1:diameter+1,1:diameter+1) = 2;

ic = 1 + radius;
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

% This step creates a matrix with the distances to the ellipse
% centroid (ic, ic, ic).  This has been optimized to minimize the
% time to generate the 3D ellipse.  The changes resulted in a
% 94% decrease in the cpu time required to generate the ellipse.

i = x0;
j = y0;
k = z0;


% this will work for first or last plane

%if ( k==ic || k==diameter + 1)
i0 = 1; i1 = diameter + 1;
j0 = 1; j1 = diameter + 1;
            a1 = [i-ic, j-ic, k-ic];
            a2 = [a^2, b^2, c^2];
            c1 = a1*M;
            c2 = c1.^2./a2;
            distance = sum(c2);
      
            
            out = distance<=  thresh; %if distance <=1 , then point in ellipse and out = 1
%end

    
    
    % If there are no pixels belonging to the ellipse on this plane,
    % then go ahead and exit out of the loop by setting k to the
    % final plane
%     if sum(sum(dist(:,:,k)<=1))==0
%         k = diameter + 1;
% end
%     % If this is not the first or last plane, then the program
%     % is smarter about which pixels it checks for ellipse pixels in
%     if k ~= ic && k ~= diameter + 1
%         d = dist(:,:,k-1) <= 1;
%         e = dist(:,:,k) <= 1;
% d1 = diff(sum(d,1)~=0);  dmin = find(d1 == 1);
% dmax = find(d1 == -1);
% if size(dmin,2) > 1; dmin = dmin(1); dmax = dmax(2); end
% e1 = diff(sum(e,1)~=0);  emin = find(e1 == 1);
% emax = find(e1 == -1);
% if size(emin,2) > 1; emin = emin(1); emax = emax(2); end
% if dmin - emin < 0;
%     j0 = emin - 1;
%     j1 = j0 + 3 + (dmax-dmin);
%     if j1 > diameter + 1; j1 = diameter + 1; end
% end
% if dmax - emax > 0;
%     j1 = emax + 1;
%     j0 = j1 - 3 - (dmax-dmin);
%     if j0 < 1; j0 = 1; end
% end
% d1 = diff(sum(d,2)~=0);  dmin = find(d1 == 1);
% dmax = find(d1 == -1);
% if size(dmin,1) > 1; dmin = dmin(1); dmax = dmax(2); end
% e1 = diff(sum(e,2)~=0);  emin = find(e1 == 1);
% emax = find(e1 == -1);
% if size(emin,1) > 1; emin = emin(1); emax = emax(2); end
% if dmin - emin < 0;
%     i0 = emin - 1;
%     i1 = i0 + 3 + (dmax-dmin);
%     if i1 > diameter + 1; i1 = diameter + 1; end
% end
% if dmax - emax > 0;
% i1 = emax + 1;
% i0 = i1 - 3 - (dmax-dmin);
% if i0 < 1; i0 = 1; end
% end
%     end
% end


% dist(i,j,k) <= 1 then (i,j,k) belongs to the ellipse