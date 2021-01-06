%% Copyright (C) 2019 José Ramom Flores das Seixas 
%%                                      <jose.ramom.flores.das.seixas@gmail.com>
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, see <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn  {} {[@var{Z}, @var{dZx}, @var{dZy}] =} zernikes_and_derivarives_cartesian_OSA (@var{x}, @var{y}, @var{n})
%% @deftypefnx {} {[@var{Z}, @var{dZx}, @var{dZy}] =} zernikes_and_derivarives_cartesian_OSA (@var{x}, @var{y}, @var{n}, @var{nan_zero})
%% Return the cartesian Zernike's pollynomials and its partial 
%% derivatives up to radial degree @var{n}, i.e. until Z[@var{n},@var{n}]
%% 
%% @var{x} is a matrix with the X coordinates of the points where the Zernike's
%% polynomials and its derivatives are computed.
%% @var{y} is a matrix with the Y coordinates of the same points.
%% @var{n} is an integer with the maximum radial degree desired.
%% @var{nan_zero} is a string that determines the values of polynomial 
%% and derived values outside the radio unit circle.
%%
%% Strictly, the polynoms are only defined for 0 <= X²+Y² <= 1.
%% If variable @var{nan_zero} = 'nan', the values of the polynomials for which
%% it is verified that (X²+Y²)>1 are set = NaN.
%% If variable @var{nan_zero} = 'zero', the values of the polynomials for which
%% it is verified that (X²+Y²)>1 are set = 0.
%%
%% @var{Z} is a 3D matrix. Each page contains a 2D matrix with the values of a 
%% Zernike's polynomial at the points defined by @var{x} and @var{y}.  
%%
%% @var{dZx} is a 3D matrix. Each page contains the values of the 
%% partial derivative in x.
%%
%% @var{dZy} is a 3D matrix. Each page contains the values of the 
%% partial derivative in y.
%% 
%% It should be noted that in standard OSA/ANSI the simple-index j starts 
%% at 0, but in octave the indices of the vectors and matrices start at 1.
%% So that page 1 of the 3D Z, dZx and dZy matrices corresponds to the 
%% single-index j = 0, and therefore to the double-index m = 0 and n = 0.
%% Page 2 corresponds to j = 1, page 3 --> j = 2, etc.
%% 
%% Example
%% @example
%% x = linspace(-1,1,101);
%% [X,Y] = meshgrid(x,x);
%% [Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (X,Y,7,'zero');
%% Z_00 = Z(:,:,1);
%%    % Z_00 is a 2D matrix with the values of Zernike's polynomial
%%    % with simple-index j = 0, and double-index m = 0 & n = 0. 
%% dZx_-24 = dZx(:,:,11);
%%    % Z_-44 is a 2D matrix with the values of the partial 
%%    % derivative in x of Zernike's polynomial with 
%%    % simple-index j = 10, and double-index m = -4 & n = 4.
%% @end example
%%
%% Run the demo to see a more complete example.
%%
%%
%% Size of @var{x} must be equal size of @var{y}.
%%
%% References:
%%
%% @enumerate
%% @item Andersen T.B., @url{https://doi.org/10.1364/OE.26.018878, "Efficient 
%%       and robust recurrence relations for the Zernike circle polynomials and
%%       their derivatives in Cartesian coordinates"}. Optic Express 26(15), 
%%       18878-18896 (2018).
%% @item Thibos, L.N, Applegate, R.A., Schwiegerling, J.T. & Webb, R., 
%%       Standards for reporting the optical aberrations of eyes. Journal of 
%%       refractive surgery, 18(5), S652-S660 (2002).
%% @end enumerate
%%
%% @seealso{zernike_osa_ansi_to_nm, zernike_cartesian, zernike_name,
%%          zernike_polar, zernike_R_poly}
%% @end deftypefn

function [Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (x, y, N, nan_zero);
  if (nargin < 3 || nargin > 4)
    print_usage ();
  elseif (~isscalar (N) || N < 1 || N ~= fix (N))
    error ("zernike_and_derivatives_cartesian_OA: Maximum radial degree must be a integer >=1");
  elseif (any (size (x) ~= size (y)))
    error ("zernike_and_derivatives_cartesian_OA: X and Y must have the same size");
  end
  r2 = x.^2 + y.^2;  out_p = find(r2>1);
  MaxJ = N*(N+3)/2;
  U1 = ones(size(x));
  switch(nan_zero)
    case "zero"
      x(out_p) = 0;
      y(out_p) = 0;
      U1(out_p) = 0;
    case "NaN"
      x(out_p) = NaN;
      y(out_p) = NaN;
      U1(out_p) = NaN;
  end
  U = zeros ([size(x) (MaxJ+1)]);    dUx = U;   dUy = U;
  Z = U;		     dZx = U;		     dZy = U;
  U(:,:,1) = U1;	dUx(:,:,1) = U1*0;	dUy(:,:,1) = U1*0;	
  Z(:,:,1) = U1;	dZx(:,:,1) = U1*0;	dZy(:,:,1) = U1*0;
  U(:,:,2) = y;	dUx(:,:,2) = U1*0;	dUy(:,:,2) = U1;
  Z(:,:,2) = 2*y;	dZx(:,:,2) = U1*0;	dZy(:,:,2) = 2*U1;
  U(:,:,3) = x;	dUx(:,:,3) = U1;	dUy(:,:,3) = U1*0;
  Z(:,:,3) = 2*x;	dZx(:,:,3) = 2*U1;	dZy(:,:,3) = U1*0;
  for k = 4:(MaxJ+1)
    [m, n] = zernike_osa_ansi_to_mn (k-1);
    switch(m)
      case -n   
	   U(:,:,k) = x.*U(:,:,(k-n))+y.*U(:,:,(k-1));
	   dUx(:,:,k) = n*U(:,:,(k-n));
	   dUy(:,:,k) = n*U(:,:,(k-1));
      case  n   
	   U(:,:,k) = x.*U(:,:,(k-(n+1)))-y.*U(:,:,(k-2*n));
	   dUx(:,:,k) = n*U(:,:,(k-(n+1)));
	   dUy(:,:,k) = -n*U(:,:,(k-2*n));
      case -1   
	   U(:,:,k) = x.*U(:,:,(k-(n+1)))+y.*U(:,:,(k-n))-y.*U(:,:,(k-(n-1)))-...
		         U(:,:,(k-2*n));
	   dUx(:,:,k) = n*U(:,:,(k-(n+1)))+dUx(:,:,(k-2*n));
	   dUy(:,:,k) = n*U(:,:,(k-n))-n*U(:,:,(k-(n-1)))+dUy(:,:,(k-2*n));
      case  1   
	   U(:,:,k) = x.*U(:,:,(k-n))+x.*U(:,:,(k-(n+1)))+y.*U(:,:,(k-(n+2)))-...
		         U(:,:,(k-2*n));
	   dUx(:,:,k) = n*U(:,:,(k-n))+n*U(:,:,(k-(n+1)))+dUx(:,:,(k-2*n));
	   dUy(:,:,k) = n*U(:,:,(k-(n+2)))+dUy(:,:,(k-2*n));
      case  0   
	   U(:,:,k) = 2*x.*U(:,:,(k-n))+2*y.*U(:,:,(k-(n+1)))-U(:,:,(k-2*n));
	   dUx(:,:,k) = 2*n*U(:,:,(k-n))+dUx(:,:,(k-2*n));
	   dUy(:,:,k) = 2*n*U(:,:,(k-(n+1)))+dUy(:,:,(k-2*n));
      otherwise 
	   U(:,:,k) = x.*U(:,:,(k-n))+y.*U(:,:,(k-(n+m+1)))+x.*U(:,:,(k-(n+1)))-...
		         y.*U(:,:,(k-(n+m)))-U(:,:,(k-2*n));
	   dUx(:,:,k) = n*U(:,:,(k-n))+n*U(:,:,(k-(n+1)))+dUx(:,:,(k-2*n));
	   dUy(:,:,k) = n*U(:,:,(k-(n+m+1)))-n*U(:,:,(k-(n+m)))+dUy(:,:,(k-2*n));
    end
    Nnm = sqrt(2*(n+1));
    if m == 0 
       Nnm = Nnm/sqrt(2);
    end
    Z(:,:,k) = Nnm * U(:,:,k);
    dZx(:,:,k) = Nnm * dUx(:,:,k);
    dZy(:,:,k) = Nnm * dUy(:,:,k);
    end
  end
%!demo
%! Nmgr = 7;                        % Maximum radial degree number, default = 7
%! Nmodos = ((Nmgr+3)*Nmgr/2)+1;    % Mode number, default = 36
%!
%! pr = 20;                         % Pupil radius, mm
%! % Let us suppose that we obtain the heights of a biconical surface
%! % for a series of points in the xy plane
%! Rx = 0.002;			          % mm⁻¹
%! px = 0.8;
%! Ry = 0.001;			          % mm⁻¹
%! py = 0.5;
%! x = pr*(2*rand(200,1) - 1);	     % mm
%! y = pr*(2*rand(200,1) - 1);	     % mm
%! surface = (Rx*x.^2 + Ry*y.^2)/(1+sqrt(1-(px*(Rx*x).^2+py*(Ry*y).^2)));
%! sag = max(max(surface));
%! surface = sag - surface;         % surface height, mm
%!
%! % Coordinates normalization
%! xnorm=x/pr;
%! ynorm=y/pr;
%! usefuldata = find(sqrt(x.^2+y.^2)<pr);     % data inside the pupil circle
%! ixnorm = xnorm(usefuldata);
%! iynorm = ynorm(usefuldata);
%! isurface = surface(usefuldata);
%!
%! % Fitting the surface			Note that ixnomr and iynorm are vectors
%! [Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (ixnorm, iynorm, Nmgr);
%! Mrec = inv(Z'*Z)*Z';             % LSQ fitting matrix
%! coefZ = Mrec*isurface;           % Fitting Zernike's coefficients
%!
%! % Reconstructed surface and partial derivatives
%! nx = linspace(-1,1,129);
%! [Xn,Yn] = meshgrid(nx,nx);
%! [Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA (Xn, Yn, Nmgr);
%! S = zeros(size(Xn));		% Note that Xn and Yn are square matrices
%! dSx = zeros(size(Xn));
%! dSy = zeros(size(Xn));
%! for m = 1:Nmodos
%!    S = S + Z(:,:,m)*coefZ(m);		% Reconstructed surface
%!    dSx = dSx + dZx(:,:,m) * coefZ(m);	% Derivative of S with respect to x
%!    dSy = dSy + dZy(:,:,m) * coefZ(m);	% Derivative of S with respect to y
%! endfor
%!
%! % Plotting the reconstructed surface and their derivatives
%! nx = linspace(-pr,pr,129);
%! [Xn,Yn] = meshgrid(nx,nx);
%! colormap ("default");
%! subplot (3,2,1:4)
%! mesh(Xn,Yn,S);
%! title("Reconstructed surface");
%! xlabel("x (mm)");      ylabel("y (mm)");        zlabel("z (mm)");
%! subplot(3,2,5)
%! imagesc(nx,nx,dSx);
%! axis square;      axis xy;       colorbar;
%! title("Derivative of Surface with respect to x");
%! subplot(3,2,6)
%! imagesc(nx,nx,dSy);
%! axis square;      axis xy;       colorbar;
%! title("Derivative of Surface with respect to y");
