%Glazed torus with rotating activity

%% Geometry
% filament is details on the slender geometry 
[filament,x] = geometry_examples(2,1,6);
filament = generate_geometry_from_xpts(x,filament);

%% Mesh
% mesh is details of surface geometry
nelm  = 100; %Number of elements
nquad = 16;
nFourier = 20; %Truncation of Fourier Series in theta
th_res = 21; % Number theta-points plotted
mesh  = generate_mesh(filament,nelm,nquad,th_res,nFourier);

%% Solutions
% Calculate SPT solutions
solution = phoretic_concentration(filament,mesh);

%Slip and swimming
solution = phoretic_slip(filament,mesh,solution);
solution = stokes_swimming(filament,mesh,solution);

% Recover solution from Fourier series expansion
[solution,mesh] = invFourierSeries(solution,filament,mesh,th_res);

%% Print swimming and rotation speeds
disp('Swimming velocity U =')
disp(solution.U)

disp('Rotational velocity Omega =')
disp(solution.Omega)