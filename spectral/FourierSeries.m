function z = FourierSeries(y,x,nF)
% Find the coefficients of a Fourier series of y up to the nF'th mode
% y is a periodic function in [0,2pi]
% cos modes from real part, sin modes from imaginary part
%%
ns = length(y(1,:)); % Number of s-points to find series for

z = zeros(nF+1,ns);

for n = 0:nF
    for k = 1:ns
        % Calculate nth Fourier mode
        cos_coeff = trapz(x,cos(n.*x).*y(:,k));
        sin_coeff = trapz(x,sin(n.*x).*y(:,k));

        % Store as cos + i sin
        z(n+1,k) = (cos_coeff + 1i*sin_coeff)./pi;
    end
end

% Zeroth order modes are half this integral
z(1,:) = z(1,:)/2;

end