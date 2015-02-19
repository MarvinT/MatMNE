function gb = gabor_fn(gamma, theta, lambd, bw, psi, sz)

% Generate a gabor patch

sigma = lambd/(pi*sqrt(log(2)/2)*(2^bw +1)/(2^bw - 1));
sigma = 5;
sigma_x = sigma;
sigma_y = sigma/gamma;


x = linspace(-sz, sz, sz);
y = linspace(sz, -sz, sz);
[xv, yv] = meshgrid(x, y);

x_theta = xv*cos(theta) +yv*sin(theta);
y_theta = -xv*sin(theta) + yv*cos(theta);
gb = exp(-0.5*((x_theta.^2/sigma_x^2) + (y_theta.^2/sigma_y^2))).*cos(2*pi/lambd .* x_theta + psi);

