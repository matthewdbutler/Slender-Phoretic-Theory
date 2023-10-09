%Planar circular arc with uniform activity

%Geometry
[filament,x] = geometry_examples(3,1,1); 
filament = generate_geometry_from_xpts(x,filament);

%Mesh
nelm  = 20; %Number of elements
nquad = 16; %Needs to be even
th_res = 21; %Doesn't affect accuracy of integrals
mesh  = generate_mesh(filament,nelm,nquad,th_res,0);

% %solution for conc at main points and offsets
solution = phoretic_concentration(filament,mesh);

%calculate slip velocity and swimming
solution = phoretic_slip(filament,mesh,solution);
solution = stokes_swimming(filament,mesh,solution);

% Recover solution from Fourier series expansion
solution = invFourierSeries(solution,filament,mesh,th_res);

%Analytical result for uniform arc
analytic = analytic_examples(3,1,mesh,filament);

%Code Validation
c0_error = 100*(analytic.c0 - solution.c0)./analytic.c0;
c1_error = 100*(analytic.c1 - solution.c1)./(analytic.c1);

%Plot
plot(mesh.seval,solution.c0)
hold on
plot(mesh.seval,analytic.c0,'--');
xlabel('s')
ylabel('conc')
legend('numeric SPT','analytic')