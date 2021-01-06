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
%% @deftypefn  {} {[@var{m}, @var{n}] =} zernike_osa_ansi_to_mn (@var{j})
%% Convert OSA/ANSI single-index @var{j} to double index @var{m} (Azimuthal 
%% degree) and @var{n} (Radial degree).
%%
%% Example
%% @example
%% [m,n] = zernike_osa_ansi_to_mn(4)
%%     @result{} [0, 2]
%% @end example
%%
%% References:
%%
%% @enumerate
%% @item Thibos, L.N, Applegate, R.A., Schwiegerling, J.T. & Webb, R., 
%%       Standards for reporting the optical aberrations of eyes. Journal of
%%       refractive surgery , 18 (5), S652-S660 (2002).
%% @item @url{https://en.wikipedia.org/wiki/Zernike_polynomials%OSA/ANSI_standard_indices,
%% OSA/ANSI standard indices of Zernike polynomials}, last retrieved on
%%       July 2109.
%% @end enumerate
%%
%% @seealso{zernike_noll_to_mn, zernikes_and_derivatives_cartesian_OSA, zernike_cartesian, zernike_name, zernike_polar, zernike_R_poly}
%%
%% @end deftypefn

function [m, n] = zernike_osa_ansi_to_mn (j)
  if (nargin ~= 1 || nargout ~= 2)
    print_usage ();
  elseif (any (j(:) < 0 | j(:) ~= fix (j(:))) || ~isvector (j))
    error ("zernike_osa_ansi_to_mn: j has to be a vector with integers >=0");
  end
  n  = ceil((-3+sqrt(9+8*j))/2);
  m  = 2*j - n .* (n+2);
end

%!test
%! [m,n]=zernike_osa_ansi_to_mn(0);
%! assert([m n],[0 0])
%! [m,n]=zernike_osa_ansi_to_mn(1);
%! assert([m n],[-1 1])
%! [m,n]=zernike_osa_ansi_to_mn(2);
%! assert([m n],[1 1])
%! [m,n]=zernike_osa_ansi_to_mn(3);
%! assert([m n],[-2 2])
%! [m,n]=zernike_osa_ansi_to_mn(4);
%! assert([m n],[0 2])
%! [m,n]=zernike_osa_ansi_to_mn(5);
%! assert([m n],[2 2])
%! [m,n]=zernike_osa_ansi_to_mn(6);
%! assert([m n],[-3 3])
%! [m,n]=zernike_osa_ansi_to_mn(7);
%! assert([m n],[-1 3])
%! [m,n]=zernike_osa_ansi_to_mn(8);
%! assert([m n],[1 3])
%! [m,n]=zernike_osa_ansi_to_mn(9);
%! assert([m n],[3 3])
%! [m,n]=zernike_osa_ansi_to_mn(10);
%! assert([m n],[-4 4])
%! [m,n]=zernike_osa_ansi_to_mn(11);
%! assert([m n],[-2 4])
%! [m,n]=zernike_osa_ansi_to_mn(12);
%! assert([m n],[0 4])
%! [m,n]=zernike_osa_ansi_to_mn(13);
%! assert([m n],[2 4])
%! [m,n]=zernike_osa_ansi_to_mn(14);
%! assert([m n],[4 4])
%! [m,n]=zernike_osa_ansi_to_mn(15);
%! assert([m n],[-5 5])
%! [m,n]=zernike_osa_ansi_to_mn(16);
%! assert([m n],[-3 5])
%! [m,n]=zernike_osa_ansi_to_mn(17);
%! assert([m n],[-1 5])
%! [m,n]=zernike_osa_ansi_to_mn(18);
%! assert([m n],[1 5])
%! [m,n]=zernike_osa_ansi_to_mn(19);
%! assert([m n],[3 5])
%! [m,n]=zernike_osa_ansi_to_mn(20);
%! assert([m n],[5 5])
%! [m,n]=zernike_osa_ansi_to_mn(21);
%! assert([m n],[-6 6])
%! [m,n]=zernike_osa_ansi_to_mn(22);
%! assert([m n],[-4 6])
%! [m,n]=zernike_osa_ansi_to_mn(23);
%! assert([m n],[-2 6])
%! [m,n]=zernike_osa_ansi_to_mn(24);
%! assert([m n],[0 6])
%! [m,n]=zernike_osa_ansi_to_mn(25);
%! assert([m n],[2 6])
%! [m,n]=zernike_osa_ansi_to_mn(26);
%! assert([m n],[4 6])
%! [m,n]=zernike_osa_ansi_to_mn(27);
%! assert([m n],[6 6])
%! [m,n]=zernike_osa_ansi_to_mn(28);
%! assert([m n],[-7 7])
%! [m,n]=zernike_osa_ansi_to_mn(29);
%! assert([m n],[-5 7])
%! [m,n]=zernike_osa_ansi_to_mn(30);
%! assert([m n],[-3 7])
%! [m,n]=zernike_osa_ansi_to_mn(31);
%! assert([m n],[-1 7])
%! [m,n]=zernike_osa_ansi_to_mn(32);
%! assert([m n],[1 7])
%! [m,n]=zernike_osa_ansi_to_mn(33);
%! assert([m n],[3 7])
%! [m,n]=zernike_osa_ansi_to_mn(34);
%! assert([m n],[5 7])
%! [m,n]=zernike_osa_ansi_to_mn(35);
%! assert([m n],[7 7])

%% vector test
%!test
%! [m,n]=zernike_osa_ansi_to_mn([2 7 15 35]);
%! assert(m,[1 -1 -5 7])
%! assert(n,[1 3 5 7])




