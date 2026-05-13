function Cl = main_kth_chw(alpha, vis)

    close all;
    if nargin < 2
        vis = false;
    end
    filename = 'airfoils\n2412.txt';
    fileid = fopen(filename);
    if fileid == -1
        error('file could not be opened');
    end
    
    coords = readmatrix(filename); %trailing edge -> along upper surface -> leading edge -> along lower -> TE
    coords = flip(coords); %trailing edge along lower surface -> le-> upper surface  -> te
    %alpha = 8; % set alpha in deg
    alpha = alpha*pi/180;
    % rot = [cos(-alpha) -sin(-alpha);
    %     sin(-alpha) cos(-alpha)]; %physically rotate the coordinates of the airfoil around its leading edge
    % coords = (rot*coords')';
    n_points = size(coords, 1);
    coords(end,:) = coords(1,:);
    
    %collocation points defined as halfway point between nodes
    coll = zeros(n_points-1,2);
    n_coll = n_points-1;
    for i = 1:(n_points-1)
        coll(i,:) = (coords(i,:)+coords(i+1,:))/2;
    end
    
    
    xb = coords(1:end-1,1); %cut off the second trailing edge point length n_coll
    %x value of nodes 2 to 1
    % xb1 = [xb(2:end); xb(1)]; %length n_coll-1
    % xb1 = coords(2:end,1);
    
    
    zb = coords(1:end-1,2);
    %z value of nodes 2 to 1
    % zb1 = [zb(2:end); zb(1)];
    % zb1 = coords(2:end,2);
    
    %angle of each panel respect to x axis
    % theta = atan2((zb1-zb), (xb1-xb));
    %length of each panel
    % s = sqrt((zb1-zb).^2 + (xb1-xb).^2);
    dx = coords(2:end,1) - coords(1:end-1,1);
    dy = coords(2:end,2) - coords(1:end-1,2);
    
    theta = atan2(dy, dx);
    s= sqrt(dx.^2 + dy.^2);
    
    xcoll=coll(:,1);
    zcoll=coll(:,2);
    LE_idx = [find(zcoll>0, 1)];
    
    Cn1 = zeros(n_coll, n_coll);
    Cn2 = zeros(n_coll, n_coll);
    
    Ct1 = zeros(n_coll, n_coll);
    Ct2 = zeros(n_coll, n_coll);
    for i=1:n_coll
        for j=1:n_coll %populating rows
            if i~=j
                Xj = xb(j);
                Zj = zb(j);
                
                %calculating coefficients
                A = -(xcoll(i)-xb(j))*cos(theta(j))-(zcoll(i)-zb(j))*sin(theta(j));
                B = (xcoll(i)-xb(j)).^2 +(zcoll(i)-zb(j)).^2;
                C = sin(theta(i)-theta(j));
                D = cos(theta(i)-theta(j));
                E = (xcoll(i)-xb(j))*sin(theta(j))-(zcoll(i)-zb(j))*cos(theta(j));
        
                F = log(1+ (s(j)^2+2*A*s(j))/B);
                G = atan2(E*s(j), B+A*s(j));
                P = ( xcoll(i)-xb(j) )*sin( theta(i)-2*theta(j) )+( zcoll(i)-zb(j) )*cos( theta(i)-2*theta(j) ) ;
                Q = ( xcoll(i)-xb(j) )*cos( theta(i)-2*theta(j) )-( zcoll(i)-zb(j) )*sin( theta(i)-2*theta(j) ) ;
                
                %log returns a complex number (even though imaginary component
                %is 0) so the real component has to be taken
                Cn2(i,j) = real(D+0.5*Q*F/s(j)-(A*C+D*E)*G/s(j)); 
                Cn1(i,j) = real(0.5*D*F+C*G-Cn2(i,j));
    
                Ct2(i,j) = real(C+0.5*P*F/s(j)+(A*D-C*E)*G/s(j));
                Ct1(i,j) = real(0.5*C*F-D*G-Ct2(i,j));
            else
                %enforcing Kutta condition
                Cn2(i,j)=1;
                Cn1(i,j)=-1;
    
                Ct2(i,j)=pi/2;
                Ct1(i,j)=pi/2;
            end
        end
    end
    
    An1 = zeros(n_coll+1, n_coll+1);
    An2 = zeros(n_coll+1, n_coll+1);
    An1(1:n_coll, 1:n_coll)= Cn1;
    An2(1:n_coll, 2:end)= Cn2;
    An = An1+An2;
    An(end,1)=1;
    An(end,end)=1;
    
    At1 =zeros(n_coll, n_coll+1);
    At2 =zeros(n_coll, n_coll+1);
    At1(1:n_coll, 1:n_coll)= Ct1;
    At2(1:n_coll,2:end)= Ct2;
    At = At1+At2;
    
    RHS = zeros(n_coll+1, 1);
    RHS(1:n_coll) = sin(theta-alpha); %uncomment to account for alpha without
    %rotating coordinates
    %RHS(1:n_coll) = sin(theta);
    gamma = An\RHS;
    
    %finding tangent velocity
    Vt = cos(theta-alpha)+At*gamma; %uncomment to account for alpha without
    %rotating coordinates
    %Vt = cos(theta)+At*gamma;
    Cp = 1-Vt.^2;
    i_positives= Cp>0;
    normals = [-sin(theta), cos(theta)];
    coll_copy = [xcoll, zcoll];
    
    
    %moved origins
    Cp_max = max(abs(Cp));
    cpbar_vecs = (Cp)/(Cp_max*8).*normals;
    coll_copy(i_positives, :) = coll_copy(i_positives, :)+cpbar_vecs(i_positives, :); %move the origins
    
    
    %finding lift coef
    %Cl = sum(-Cp.*s.*cos(theta));
    Cl = 2 * sum(gamma(1:end-1) .* s);
    
    if vis
        figure()
        hold on
        plot(xcoll, zcoll, 'or')
        plot(xb,zb, '-ob')
        %plot(cpoutline_x,cpoutline_z)
        quiver(coll_copy(:,1), coll_copy(:,2), -cpbar_vecs(:,1), -cpbar_vecs(:,2), 0, 'MaxHeadSize', 0.2)
        axis equal
        
        figure()
        hold on
        plot(xcoll(1:LE_idx), Cp(1:LE_idx), '-.r', DisplayName='lower surface');
        plot(xcoll(LE_idx:end), Cp(LE_idx:end), '-.b', DisplayName='upper surface');
        set(gca, 'YDir', "reverse", ...
            'Ylabel', ylabel('C_p'), ...
            'Xlabel', xlabel('x/c'));
        legend();
    end

end