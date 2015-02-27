## Copyright (C) 2002 André Carezia <acarezia@uol.com.br>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {} chebwin (@var{m})
## @deftypefnx {Function File} {} chebwin (@var{m}, @var{at})
##
## Return the filter coefficients of a Dolph-Chebyshev window of length @var{m}.
## The Fourier transform of the window has a stop-band attenuation of @var{at}
## dB.  The default attenuation value is 100 dB.
##
## For the definition of the Chebyshev window, see
##
## * Peter Lynch, "The Dolph-Chebyshev Window: A Simple Optimal Filter",
##   Monthly Weather Review, Vol. 125, pp. 655-660, April 1997.
##   (http://www.maths.tcd.ie/~plynch/Publications/Dolph.pdf)
##
## * C. Dolph, "A current distribution for broadside arrays which
##   optimizes the relationship between beam width and side-lobe level",
##   Proc. IEEE, 34, pp. 335-348.
##
## The window is described in frequency domain by the expression:
##
## @example
## @group
##          Cheb(m-1, beta * cos(pi * k/m))
##   W(k) = -------------------------------
##                 Cheb(m-1, beta)
## @end group
## @end example
##
## with
##
## @example
## @group
##   beta = cosh(1/(m-1) * acosh(10^(at/20))
## @end group
## @end example
##
## and Cheb(m,x) denoting the m-th order Chebyshev polynomial calculated
## at the point x.
##
## Note that the denominator in W(k) above is not computed, and after
## the inverse Fourier transform the window is scaled by making its
## maximum value unitary.
##
## @seealso{kaiser}
## @end deftypefn

function w = chebwin (m, at)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (! (isscalar (m) && (m == fix (m)) && (m > 0)))
    error ("chebwin: M must be a positive integer");
  elseif (nargin == 1)
    at = 100;
  elseif (! (isscalar (at) && isreal (at)))
    error ("chebwin: AT must be a real scalar");
  endif

  if (m == 1)
    w = 1;
  else
    ## beta calculation
    gamma = 10^(-at/20);
    beta = cosh(1/(m-1) * acosh(1/gamma));
    ## freq. scale
    k = (0:m-1);
    x = beta*cos(pi*k/m);
    ## Chebyshev window (freq. domain)
    p = cheb(m-1, x);
    ## inverse Fourier transform
    if (rem(m,2))
      w = real(fft(p));
      M = (m+1)/2;
      w = w(1:M)/w(1);
      w = [w(M:-1:2) w]';
    else
      ## half-sample delay (even order)
      p = p.*exp(j*pi/m * (0:m-1));
      w = real(fft(p));
      M = m/2+1;
      w = w/w(2);
      w = [w(M:-1:2) w(2:M)]';
    endif
  endif

  w = w ./ max (w (:));

endfunction

%!assert (chebwin (1), 1)
%!assert (chebwin (2), ones (2, 1))

%% Test input validation
%!error chebwin ()
%!error chebwin (0.5)
%!error chebwin (-1)
%!error chebwin (ones (1, 4))
%!error chebwin (1, 2, 3)
