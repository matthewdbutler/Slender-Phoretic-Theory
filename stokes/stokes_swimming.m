function solution = stokes_swimming(filament,mesh,solution)
% Stokes solver for swimming given the slip velocities

% Aim is to invert Ax=b
% A is (3N+6)x(3N+6) matrix
% x is (3N+6)x1 vector of N forces/length (x 8pi), plus the swimming 
% translational and angular velocities
% b is (3N+6)x1 vector of N theta-averaged slip velocities and 6 zeros

% The 3Nx3N part of A is composed of 3x3 submatrices Mij which are:
% Int(G-J) over segment - Int(J) over everything else + L , when i=j
% Integral over sdum of jth segment of G(s_i,sdum), otherwise

% G = (I + R0_hat R0_hat)/|R0|
% J = (I + t t)/|sdum-s|
% L = log(4(1-s^2)/(eps^2 rho^2)) (I + t t) + I - 3 t t
% The log term does not have (1-s^2) when shape is looped, and q=sdum-s

% From Taylor expansion, find integral over segment
% Int(G-J) = K^2 (nn/4 - 3tt/8 - I/8) l + O(l^2)
% Integrating away from there, also have
% Int(G) = log(4(1-s^2)/l^2) (I + tt)

neval = length(mesh.seval); % number of s-points
nquad = length(mesh.squad(:,1)); % number of quatrature points

xind  =         1:nquad;
yind  =   nquad+1:2*nquad;
zind  = 2*nquad+1:3*nquad;

A = zeros(3*neval+6);
for i=1:neval
    
    % 3x3 tt matrix
    tangent_matrix = mesh.tangent_eval(:,i)*mesh.tangent_eval(:,i)';
    
    % Calculate different log terms
    if filament.isopen == 0 % closed filament

        log_bit = log(4./filament.epsilon^2./mesh.rho_eval(i).^2);
        
        q1 = mesh.squad - mesh.seval(i);
        q  = q1;
        q(q1 > 1)  = q(q1 > 1)  - 2;
        q(q1 < -1) = q(q1 < -1) + 2;

    else % open filament
        log_bit = log(4*(1-mesh.seval(i).^2)./filament.epsilon^2./mesh.rho_eval(i).^2);
        
        %|s - s'| denominator for regularisation
        q = mesh.squad - mesh.seval(i);
    end
    
    %Quadrature vectors
    dx = mesh.xquad(xind,:) - mesh.xeval(1,i);
    dy = mesh.xquad(yind,:) - mesh.xeval(2,i);
    dz = mesh.xquad(zind,:) - mesh.xeval(3,i);
    
    %R_0 in equation (distance in integrand)
    r2 = dx.^2 + dy.^2 + dz.^2;
    r = sqrt(r2);
    
    for j=1:neval
        
        if i==j % M = int(G-J) - int(J) + L
            L = log_bit.*(eye(3) + tangent_matrix) + eye(3) - 3.*tangent_matrix;
            
            % Make 3x3 matrices for each quadrature point and sum
            R0_quad = [dx(:,j)';dy(:,j)';dz(:,j)'];
            GJintegrand = zeros(3,3,nquad);
            for quad = 1:nquad
                GJintegrand(:,:,quad) = eye(3)./r(quad,j) + R0_quad(:,quad)*R0_quad(:,quad)'./r(quad,j).^3;
                GJintegrand(:,:,quad) = GJintegrand(:,:,quad) - (eye(3)+tangent_matrix)./abs(q(quad,j));
                GJintegrand(:,:,quad) = GJintegrand(:,:,quad).*mesh.quad_weights(quad,j);
            end
            GJinner = sum(GJintegrand,3);
            
            % integrate J over all segments except i  
            Jmagnitude = 0;
            for k=1:neval
                if k~=i   
                    Jintegrand = mesh.quad_weights(:,k)./abs(q(:,k));
                    Jmagnitude = Jmagnitude + sum(Jintegrand);
                end
            end
            Jouter = Jmagnitude.*(eye(3)+tangent_matrix);
            
            M = GJinner - Jouter + L;
            
        else % M = int(G)
            
            % Make 3x3 matrices for each quadrature point and sum
            R0_quad = [dx(:,j)';dy(:,j)';dz(:,j)'];
            Gintegrand = zeros(3,3,nquad);
            for quad=1:nquad
                Gintegrand(:,:,quad) = eye(3)./r(quad,j) + R0_quad(:,quad)*R0_quad(:,quad)'./r(quad,j).^3;
                Gintegrand(:,:,quad) = Gintegrand(:,:,quad).*mesh.quad_weights(quad,j);
            end
            
            M = sum(Gintegrand,3);

        end
        
        A(3*i-2:3*i,3*j-2:3*j) = M;
  
    end
    
    % Add identity sub-matrices
    A(3*i-2:3*i,3*neval+1:3*neval+3) = -eye(3);
    A(3*neval+1:3*neval+3,3*i-2:3*i) = eye(3);
    
    % Add Levi-Civita sub-matricesjanus
    eps_dot_r = -[cross([1;0;0],mesh.xeval(:,i)) cross([0;1;0],mesh.xeval(:,i)) cross([0;0;1],mesh.xeval(:,i))];
    A(3*i-2:3*i,3*neval+4:3*neval+6) = eps_dot_r;
    A(3*neval+4:3*neval+6,3*i-2:3*i) = eps_dot_r;
end

vslip0 = [solution.vlead_x(1,:); solution.vlead_y(1,:); solution.vlead_z(1,:)];

% Make vector
b = [reshape(vslip0,[],1); zeros(6,1)];

% If filament is straight, remove Omega_x contributions (singular solution)
if max(abs(mesh.curvature_eval))<1e-6
    A(:,3*neval+4) = [];
    
    % Invert to solve
    x = A\b;
    
    solution.U = x(3*neval+1:3*neval+3);
    solution.Omega = [0;x(3*neval+4:3*neval+5)];
else
    % Invert to solve
    x = A\b;
    
    solution.U = x(3*neval+1:3*neval+3);
    solution.Omega = x(3*neval+4:3*neval+6); 
end

end
