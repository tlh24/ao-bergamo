Nmgr = 7;                        % Maximum radial degree number, default = 7
Nmodos = ((Nmgr+3)*Nmgr/2)+1;    % Mode number, default = 36

pr = 20;                         % Pupil radius, mm
% Let us suppose that we obtain the heights of a biconical surface
% for a series of points in the xy plane
Rx = 0.002;			          % mm⁻¹
px = 0.8;
Ry = 0.001;			          % mm⁻¹
py = 0.5;
x = pr*(2*rand(200,1) - 1);	     % mm
y = pr*(2*rand(200,1) - 1);	     % mm
surface = (Rx*x.^2 + Ry*y.^2)./(1+sqrt(1-(px.*(Rx*x).^2+py.*(Ry*y).^2)));
sag = max(max(surface));
surface = sag - surface;         % surface height, mm

% Coordinates normalization
xnorm=x/pr;
ynorm=y/pr;
usefuldata = find(sqrt(x.^2+y.^2)<pr);     % data inside the pupil circle
ixnorm = xnorm(usefuldata);
iynorm = ynorm(usefuldata);
isurface = surface(usefuldata);

% Fitting the surface			Note that ixnomr and iynorm are vectors
[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (ixnorm, iynorm, Nmgr, "NaN");
% Mrec = inv(Z'*Z)*Z';             % LSQ fitting matrix
% want to fit Z * coefZ = isurface
Zs = squeeze(Z); 
coefZ = Zs \ isurface;
% coefZ = Mrec*isurface;           % Fitting Zernike's coefficients

% Reconstructed surface and partial derivatives
nx = linspace(-1,1,129);
[Xn,Yn] = meshgrid(nx,nx);
[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (Xn, Yn, Nmgr, "NaN");
S = zeros(size(Xn));		% Note that Xn and Yn are square matrices
dSx = zeros(size(Xn));
dSy = zeros(size(Xn));
for m = 1:Nmodos
	S = S + Z(:,:,m)*coefZ(m);		% Reconstructed surface
	dSx = dSx + dZx(:,:,m) * coefZ(m);	% Derivative of S with respect to x
	dSy = dSy + dZy(:,:,m) * coefZ(m);	% Derivative of S with respect to y
end

% Plotting the reconstructed surface and their derivatives
nx = linspace(-pr,pr,129);
[Xn,Yn] = meshgrid(nx,nx);
colormap ("default");
subplot (3,2,1:4)
mesh(Xn,Yn,S);
title("Reconstructed surface");
xlabel("x (mm)");      ylabel("y (mm)");        zlabel("z (mm)");
subplot(3,2,5)
imagesc(nx,nx,dSx);
axis square;      axis xy;       colorbar;
title("Derivative of Surface with respect to x");
subplot(3,2,6)
imagesc(nx,nx,dSy);
axis square;      axis xy;       colorbar;
title("Derivative of Surface with respect to y");
