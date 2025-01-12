function [C1,C2,C3] = find_multidim_contour(mu,Sigma,CI_Z)
% Given mean vector (mu) and covariance matrix (Sigma), find the confidence
% intervals (CIs) for 90%, 95%, and 99%.
% Output: C is a nx-by-num array, where nx is the dimension of state
dim = size(mu,2);
if dim==2
    elpt = ellipsedata(Sigma,mu,100,CI_Z);
    C1(1,:) = elpt(:,1)';
    C1(2,:) = elpt(:,2)';
    C2(1,:) = elpt(:,3)';
    C2(2,:) = elpt(:,4)';
    C3(1,:) = elpt(:,5)';
    C3(2,:) = elpt(:,6)';
elseif dim==4
    C1 = contour4D(CI_Z(1),mu,Sigma);
    C2 = contour4D(CI_Z(2),mu,Sigma);
    C3 = contour4D(CI_Z(3),mu,Sigma);
else
    error('Function `find_multidim_contour` can only be used for 2-D or 4-D cases')
end

function x = contour4D(c,mu,Sigma)
    [PD,PV] = eig(Sigma);
    PV = diag(PV);
    for i = 1:4
        a(i) = c*sqrt(PV(i));
    end
    num_pts = 30;
    t1 = linspace(0,pi,num_pts);
    t2 = linspace(0,2*pi,num_pts);
    t3 = linspace(0,2*pi,num_pts);
    counter = 1;
    % Construct ellipse
    for i1 = 1:num_pts
        for i2 = 1:num_pts
            for i3 = 1:num_pts
                x(1,counter) = a(1)*sin(t1(i1)).*sin(t2(i2)).*cos(t3(i3));
                x(2,counter) = a(2)*sin(t1(i1)).*sin(t2(i2)).*sin(t3(i3));
                x(3,counter) = a(3)*sin(t1(i1)).*cos(t2(i2));
                x(4,counter) = a(4)*cos(t1(i1));
                counter = counter+1;
            end
        end
    end
    x = x'*PD';
    x = x + repmat(mu,size(x,1),1);
    x = x';
end
end