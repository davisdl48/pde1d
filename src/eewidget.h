/*
    Euler Explicit Method - First Order
    Copyright (C) 2012  Davis Family <email>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#ifndef EEWIDGET_H
#define EEWIDGET_H

#include "solvwidget.h"

/// First Order Euler Explicit
/**
The Euler Explicit method is one the simplest numerical initial value Ordinary 
Differential Equation solvers.  It is based on an approximation to the derivative
where the change in independent variable is small but not approaching zero.

The derivative or slope of a function \f$ u \f$ at a location,\f$ x_i \f$ , and time,
\f$ t_n \f$ , with respect to time is defined as 
\f[ \frac{ \partial u }{ \partial t } = \lim_{\delta t \to 0 } 
\frac { u(x_i,t_n + \delta t ) - u(x_i , t_n ) }{\delta t} \f]
where \f$ x_i \f$ is the i<sup>th</sup> x point and \f$ t_n \f$ is the n<sup>th</sup> time step.  By replacing 
\f$ \delta t \f$ with a small \f$ \Delta t \f$ on the right side and solving for 
\f$ u(x_i,t_n + \Delta t ) \f$ an approximation for the value of \f$ u \f$ at the new time
can be obtained.
\f[ u(x_i,t_n + \Delta t ) = u(x_i,t_{n}) + \Delta t \frac{ \partial u }{ \partial t } \f]

\par Time Component 

The Euler Explicit method uses the partial derivative with respect to time at the current
time to predict the value at the next time step, that is
\f[ u_i^{n+1} = u_i^{n} + \Delta t u_t \f]
where \f$ u_i^{n+1} \f$ is notation for \f$ u(x_i,t_n + \Delta t) \f$ and \f$ u_t \f$
is shorthand notation for \f$ \frac{ \partial u }{ \partial t } \f$.  Thoughout this discusion
subscripts of t or x will represent partial derivatives while supscripts of i or k will represent 
indecies.
 
The temporal partial derivative, \f$ u_t \f$ , is solved from the partial differential equation(PDE).
Where the spacial derivatives are also approximated from finite differences.
 
\par Spacial Components
 
The spacial components used here all assume uniform spacing, \f$ \Delta x \f$.
The first spacial derivatives \f$ \frac{\partial u}{\partial x} \equiv u_x \f$ are
approximated as upwind differences based on the sign of it's coefficient.  That is
\f$ c u_x \approx c \left( \frac { u_i - u_{i-1} }{\Delta x} \right) \f$ if c is positive and
\f$ c u_x \approx c \left( \frac { u_{i+1} - u_{i} }{\Delta x} \right) \f$ if c is negative.
Another way to express this is 
\f[ c u_x \approx \frac{c}{2 \Delta x} ( u_{i+1} - u_{i-1} ) + \frac{|c|}{2 \Delta x} 
 ( -u_{i+1} + 2 u_i -u_{i-1} ) \f] 
where \f$ |c| \f$ is the absolute value of \f$ c \f$.
 
For second derivative terms, a central difference is used, 
 \f$ \nu u_{xx} \approx \nu \frac { u_{i+1} - 2 u_i + u_{i-1} }{ \Delta x^2 } \f$.
   
\par Advection/Diffusion Equation
 
Substituting the above approximations into the Advection/Diffusion equation, 
\f$ u_t + c u_x = \nu u_xx \f$ , results in
\f{eqnarray*}{ u_i^{n+1} = u_i^{n} + \Delta t &(& -\frac{c}{2 \Delta x} ( u_{i+1} - u_{i-1} ) \\
&& \quad - \frac{|c|}{2 \Delta x}  ( -u_{i+1} + 2 u_i -u_{i-1} ) 
 + \nu \frac { u_{i+1} - 2 u_i + u_{i-1} }{ \Delta x^2 } \\
   &)& \f}
or combining a few terms give the difference equation
\f{eqnarray*}{   u_i^{n+1} = u_i^{n} + \frac{\Delta t}{\Delta x}  &(& \frac{-c}{2} ( u_{i+1} - u_{i-1} ) \\
 && \quad +( \frac{\nu}{\Delta x} + \frac{|c|}{2} )  ( u_{i+1} - 2 u_i + u_{i-1} ) \\
   &)& \f}.  
\note The term with \f$ |c| \f$ is added to the diffusion term, \f$ \nu \f$ , which shows that 
upwinding has the same form as adding additional diffusion.  The 'nulimit' option described
below reduces \f$ \nu \f$ to account for this dissipation.
 
\par Finite Volume Form

In one dimension we will call the line from midpoints \f$ x_{i-1/2} \f$ to \f$ x_{i+1/2} \f$ 
the volume.   Form a box in space and time with corners at \f$ (x_{i-1/2},t_n)
, (x_{i+1/2},t_n), (x_{i-1/2},t_{n+1}), (x_{i+1/2},t_{n+1}) \f$.  Then rewritting the 
advection diffusion equation as \f$ u_t + f_x = 0 \f$ where \f$ f = c u - \nu u_x \f$, since 
\f$ \nu \f$ is a constant.  Integrating the above equation in time gives 
\f$ u(x,t_n + \Delta t) - u(x,t_n) = - \int_{t_n}^{t_n + \Delta t}f_x dt \f$.  Then the
'volume' average of \f$ u \f$ , \f$ \bar u \f$ , at the new time level becomes
\f{eqnarray*}{ 
   \overline{u^{n+1}} &=& \frac{1}{\Delta x} \int_{x_{i-1/2}}^{x_{i+1/2}} u(x,t_n + \Delta t)\, dx \\
 &=& \frac {1}{\Delta x} \int_{x_{i-1/2}}^{x_{i+1/2}} u(x,t_n) - \int_{t_n}^{t_n + \Delta t}f_x \,dt \, dx \\
 &=& \frac {1}{\Delta x} \int_{x_{i-1/2}}^{x_{i+1/2}} u(x,t_n)\,dx 
 - \int_{t_n}^{t_n + \Delta t}\int_{x_{i-1/2}}^{x_{i+1/2}}f_x \,dx \, dt
\f}
The integrals can be reversed because we assume that \f$ f_x \f$ is continuous
and differentiable.  Also, from the Fundamental theorem of calculus (the divergence theorem )
 \f[ \int_{x_{i-1/2}}^{x_{i+1/2}} f_x \, dx = f(x_{i+1/2},t) -f(x_{i-1/2},t) \f] which can
be substituted into the above equation to yield
 \f[ \overline{u_i^{n+1}} = \overline{ u_i^n} - \frac {1}{\Delta x}
 \int_{t_n}^{t_n + \Delta t} f(x_{i+1/2},t) -f(x_{i-1/2},t) dt \f]
 
Then assuming \f$ f \f$ is constant with respect to time over the interval from \f$ t \f$ to 
\f$ t + \Delta t \f$ and solving for \f$ \bar{u} \f$ we obtain
\f{eqnarray*}{
    \overline{u_i^{n+1}} &=& \overline{u_i^n} + \frac{ \Delta t }{\Delta x} \left( f(x_{i-1/2}) - f(x_{i+1/2}) \right) \\
    &=& \overline{u_i^n} + \frac{ \Delta t }{\Delta x} \left( c u_{i-1/2} - \nu u_x|_{i-1/2} 
    -  c u_{i+1/2} +  \nu u_x|_{i+1/2}  \right)
\f}
  
Again upwinding for the advection term gives 
\f[ c u_{i-1/2} = \frac{c}{2} ( u_{i-1} + u_i ) + \frac{|c|}{2} ( u_{i-1} - u_i ) \f] and
\f[ c u_{i+1/2} = \frac{c}{2} ( u_{i} + u_{i+1} ) + \frac{|c|}{2} ( u_i - u_{i+1} ) \f].
Also, use the approximation \f[ u_x|_{i-1/2} = \frac{u_i - u_{i-1}}{\Delta x} \f] and
\f[ u_x|_{i+1/2} = \frac{u_{i+1} - u_i }{\Delta x} \f].  Dropping the overbar for the
average of \f$ u \f$  gives
\f[ u_i^{n+1} =  u_i^n + \frac{\Delta t }{\Delta x} \left( \frac{c}{2} ( u_{i-1} - u_{i+1} ) \right)   + \frac{\Delta t }{\Delta x^2} 
  \left( \nu + \frac{|c| \Delta x}{2} \right) (u_{i+1} - 2 u_i + u_{i-1}) \f].
This is the same result as the finite difference equation above.

\note For conservation law PDEs the derivation usually involves the reverse process.  That is 
the finite volume form of the conservation equations is converted into the finite difference form
using the divergence theorm.  The finite volume form is then a more natural expression of the conservation
equations.  However, accuracy and stability analysis is simpler in the finite difference form. 
 
Accuracy ( Convergence )
------------------------
 
To verify the accuracy of this equation substitute the Taylor Series expansion of \f$ u \f$
about \f$ x_i \f$.
 \f[ u_{i-1} = u(x_i - \Delta x,t_n) = u(x_i) - \Delta x \, u_x(x_i,t_n)
 +  \frac {\Delta x^2}{2} u_{xx}(x_i,t_n) - \frac {\Delta x^3}{3!} u_{xxx}(x_i,t_n) + R_k(x) \f]
 \f[ u_{i+1} = u(x_i + \Delta x,t_n) = u(x_i) + \Delta x u_x(x_i,t_n)
 +  \frac {\Delta x^2}{2} u_{xx}(x_i,t_n) + \frac {\Delta x^3}{3!} u_{xxx}(x_i,t_n) + R_k(x) \f]
where \f$ R_k(x)\f$ is a remander or error term and for functions \f$ u \f$ that are
 \f$ k+1 \f$ times differentiable on the open interval and continuous on the closed interval between
 \f$ x\f$ and \f$ x \pm \Delta x  \f$ is
 \f[ R_k (x) = \frac {u_{x^(k+1)}(\xi_L)}{(k+1)!} \Delta x^{k+1} \f]
for some \f$ \xi_L \f$ in the interval between \f$ x \f$  and \f$ x \pm \Delta x \f$.  In the
above equations \f$ k = 3 \f$ and \f$ u_{x^(k+1)} \f$ is the (k+1)<sup>th</sup> partial derivative of 
\f$ u \f$ with respect to \f$ x \f$.
 
With these relations, the difference terms become
\f[ u_{i+1} - u_{i-1} = 2 \Delta x u_x + 2 \frac { \Delta x^3}{3!} u_{xxx} + R_k(x) \f]
and
\f[ u_{i+1} - 2 u_i + u_{i-1} = \Delta x^2 u_{xx} + 2 \frac {\Delta x^4}{4!} u_{x^{(4)}} + R_k(x) \f]
and for the time derivative
\f[ u_i^{n+1} =   u_i^n + \Delta t u_t + \frac {\Delta t^2}{2} u_{tt}
 + \frac {\Delta t^3}{3!} u_{ttt} + R_k(t) \f].
Then substituting into the differnce equation above gives
\f{eqnarray*}{  u_i^n + \Delta t u_t &+& \frac {\Delta t^2}{2} u_{tt} + \frac {\Delta t^3}{3!} u_{ttt} \\
 &=& u_i^{n} + \frac{\Delta t}{\Delta x}  
    \left(    -\frac{c}{2} ( 2 \Delta x u_x + 2 \frac { \Delta x^3}{3!} u_{xxx} ) 
 + \left( \frac{\nu}{\Delta x} + \frac{|c|}{2} \right)  
 ( \Delta x^2 u_{xx} + 2 \frac {\Delta x^4}{4!} u_{x^{(4)}} )
    \right) + ... \f} .
Collecting some terms, keeping only terms up to first order deltas gives
    \f[  u_t + c u_x - \nu u_{xx}  = \frac{|c| \Delta x}{2} u_{xx} - \frac {\Delta t}{2} u_{tt} \f]
which approaches the advection/diffusion equation as \f$ \Delta x \f$ and \f$ \Delta t \f$ approach zero.
Because the error term on the right hand side involves first order deltas the method is first order
in time and space.
 
The 'nulimit' option reduces the viscosity by the aparent viscosity do to upwinding through
the absolute value of wave speed, \f$ |c| \f$.
 \f[ \nu_{reduced} = max( 0.0, \nu - \frac{|c| \Delta x}{2} ) \f]
    
Stability
---------

Stability is one of the major cornerstones to numerical methods.  Indeed if a method is not 
stable it will not give the correct solution to the partial differential equation.  When a 
numerical method gives the correct solution is called consistancy.  The are two requirements
for consistancy, convergence and stability.  Convergence for the Explicit Euler method was 
shown above under accuracy.  Convergence indicates that the numerical difference equation
approaches the partial differential equation as the step size in both time and space approach
zero.  However, convergence may be limitted by roundoff or truncation errors introduced by the 
finite precision of computer hardware.

Stability analysis for non-linear equation and including boundary conditions becomes very 
difficult.  However, for linear systems with periodic bounday conditions the stability analysis
can give insight into the stability requirements of the numerical method.  Linear stability 
is a requirement for stability but it does not in the general case, prove stability.  Non-linear 
equations can be linearized for stability analysis but non-linear interactions may still be unstable
which would not be indicated in the linear analysis.

For linear partial differential systems the solution is a sum of solutions that depend on the
initial and boundary conditions.  The solutions can be based on sum of sine and cosine waves 
since derivatives of sines and cosines are sines and cosines.  Again, since the system is linear
we can examine one component, a sine and a cosine wave at a frequency, independently.  If any
component is unstable, the method is unstable.

\par Linear Diffusion Stability

First, examine the diffusion part of the equation
\f{eqnarray*}{
  u_t &=& \nu \: u_{xx} \\[3mm]
  \frac{ u_i^{n+1} - u_i^n }{\Delta t} &=& \nu \: \frac{ u_{i-1}^n - 2 u_i^n + u_{n+1}^n }{\Delta x^2}
\f}
where one component of \f$ {u_k^n} = a \sin (\omega_k x) + b \cos (\omega_k x) \f$.  Where \f$ k \f$
indicates a solution component that is continuous in \f$ x \f$ at a frequency of \f$ \omega_k = k \, \omega \f$.
Then using the trigonometric angle addition formulas for sine and cosine

\see http://mathworld.wolfram.com/TrigonometricAdditionFormulas.html , http://www.themathpage.com/atrig/sum-proof.htm
or http://en.wikipedia.org/wiki/Proofs_of_trigonometric_identities for angle sum formulas and proofs

\f{eqnarray*}{
  u_k^n |_{i+1} &=& a \sin (\omega_k x + \omega_k \, \Delta x) + b \cos (\omega_k x + \omega_k\, \Delta x) \\
  &=& a ( \cos (\omega_k x ) \sin (\omega_k \, \Delta x) + \sin (\omega_k x) \cos (\omega_k \, \Delta x ) ) \\
  && \   + b ( \cos (\omega_k x ) \cos (\omega_k \, \Delta x) - \sin (\omega_k x) \sin (\omega_k \, \Delta x ) ) \\[3mm]
  u_k^n |_{i-1} &=& a ( - \cos (\omega_k x ) \sin (\omega_k \Delta x) + \sin (\omega_k x) \cos (\omega_k  \Delta x ) ) \\
  && \  + b ( \cos (\omega_k x ) \cos (\omega_k \, \Delta x) + \sin (\omega_k x) \sin (\omega_k \, \Delta x ) ) 
\f}
results in
\f{eqnarray*}{
u_{k,i-1}^n - 2 u_{k,i}^n + u_{k,i+1}^n &=& 2 a \; \sin (\omega_k \, x) \cos ( \omega_k \, \Delta x ) +
2 b \, \cos (\omega_k \, x) \cos ( \omega_k \, \Delta x )\\
&& \ -2 a \, \sin (\omega_k \, x) - 2 b \, \cos (\omega_k \, x) \\
&=& 2 (  a \, \sin ( \omega_k \, x) ( \cos (\omega_k \, \Delta x) - 1 ) \\
&& \ + b \, \cos ( \omega_k \, x) ( \cos (\omega_k \, \Delta x) - 1 ) \\
&=& 2 u_k^n ( \cos (\omega_k \, \Delta x) - 1 )
\f}
Substituting into the PDE gives
\f{eqnarray*}{ {u_k^{n+1}} &=& {u_k^n} + 
2 \frac{ \nu \, \Delta t}{ \Delta x^2} {u_k^n} ( \cos ( \omega_k \, \Delta x ) -1 )\\
&=& {u_k^n} \left( 1 + 2 \frac{ \nu \, \Delta t}{ \Delta x^2}
( \cos ( \omega_k \, \Delta x ) -1 ) \right)
\f}
This indicates the component is multiplied by 
\f$ g = 1 + 2 \frac{ \nu \, \Delta t}{ \Delta x^2}( \cos ( \omega_k \, \Delta x ) -1 ) \f$
each time step.  If \f$ |g| > 1 \f$ the magnitude of that component will
increase without bound as <em>n</em> approaches infinity.
Since \f$ -2 \le ( \cos ( \omega_k \, \Delta x ) -1 ) \le 0 \f$ ,
\f[ 1 - 4 \frac{ \nu \, \Delta t}{ \Delta x^2} \le g \le 1 \; .\f]
Then \f$ |g| \le 1 \f$ requires
\f[ 1 - 4 \frac{ \nu \, \Delta t}{ \Delta x^2} \ge -1 \f]
or, since all terms, \f$ ( \nu , \Delta t , \Delta x ) \f$, are positive,
\f$ \frac{ \nu \, \Delta t}{ \Delta x^2} \le \frac{1}{2} \f$ or
\f$ \Delta t \le \frac{\Delta x^2}{ \nu \, 2} \f$.  This indicates a conditional stability
where the time step is limitted by stability constraints.

\note A series of values at equally spaced points can be interpolated with a series of sine and
cosine functions.  The process of converting from point values to sine and cosine functions is
a Discrete Fourier Transform.  A special form of Discrete Fourier Transform is the Fast Fourier Transform
that is both computationally efficient and less susceptible to round off errors than straight forward 
Discrete Fourier Transform computation.

That is for each \f$ i \f$  
\f[ u(x_i) = \sum_{k=0}^K { a_k \, \sin ( \omega_k \, x_i ) + b_k \, \cos ( \omega_k \, x_i )} \f]
or introducing complex numbers with \f$ I = \sqrt{ -1} \f$, and Euler's Formula,
\f$ e^{I \, x} = \cos(x) + I \, \sin(x) \f$
\f{eqnarray*}{
  u(x_i) &=& \sum_{k=-K}^K c_k e^{\omega_k \, x_i} \\
&=& \sum_{k=-K}^K {c_k ( \, I \, \sin ( \omega_k \, x_i ) +  \, \cos ( \omega_k \, x_i ) ) }
\f}
where \f$ c_k = \frac{b_k}{2} - I \, \frac{a_k}{2} ,\  a_{-k} = -a_k \f$  and \f$ b_{-k} = b_k \f$.

To show that the above equations for \f$ u(x_i) \f$ are equivalent look at the sum of
a \f$ \pm k \f$ pair.
\f{eqnarray*}{
u_k + u_{-k} &=& c_k ( \, I \, \sin ( \omega_k \, x_i ) +  \, \cos ( \omega_k \, x_i ) ) \\
&& \qquad  + c_{-k} ( \, I \, \sin ( -\omega_k \, x_i ) +  \, \cos ( - \omega_k \, x_i ) ) \\
&=& (\frac{b_k}{2} - I \, \frac{a_k}{2} ) ( \, I \, \sin ( \omega_k \, x_i ) +  \, \cos ( \omega_k \, x_i ) ) \\
&& \qquad  + (\frac{b_k}{2} + I \, \frac{a_k}{2} ) ( \, -I \, \sin ( \omega_k \, x_i ) +  \, \cos ( \omega_k \, x_i ) ) \\
&=& b_k \, \cos ( \omega_k \, x_i ) - I^2 a_k \sin ( \omega_k \, x_i ) + 0 I
\f}
which, with the addition of \f$ c_0 = b_0 \f$, results in the first equation.

\par Advection Stability
 
\par Burgers' Equation

Burgers' equation, \f$ u_t + \frac{1}{2} ( u^2 )_x = \nu u_{xx} \f$ or 
\f$ u_t +  u \, u_x = \nu u_{xx} \f$ is a non-linear PDE.  The wave speed is 
a function of the dependent variable, \f$ u \f$.  The greater \f$ u \f$ the faster the movement.
An initial line with a negative slope ( u decreasing with x ) will steepen as time moves forward
while an initial line with a positive slope will flatten as time moves forward.
  
With the finite difference option, the upwinding is implimented as

\f$ u_i^{n+1} = u_i^n + \Delta t \left( \frac{ u_{i-1}^2 - u_i^2}{ 2 \Delta x} 
+ \nu \frac {u_{i-1} - 2 u_i + u_{i+1}}{\Delta x^2} \right) \f$ - if \f$ u_i > 0 \f$ ,

\f$ u_i^{n+1} = u_i^n + \Delta t \left( \frac{ u_{i}^2 - u_{i+1}^2}{ 2 \Delta x} 
+ \nu \frac {u_{i-1} - 2 u_i + u_{i+1}}{\Delta x^2} \right) \f$ - if \f$ u_i < 0 \f$ and

\f$ u_i^{n+1} = u_i^n + \Delta t \left(  
+ \nu \frac {u_{i-1} - 2 u_i + u_{i+1}}{\Delta x^2} \right) \f$ - if \f$ u_i = 0 \f$ .

Applying the Taylor Series expansion to all terms for the case \f$ u_i > 0 \f$ and rearranging gives
\f[ u_t + u u_x - \nu u_{xx} = \frac { \Delta x}{2} ( (u_x)^2 + u u_{xx}) - \frac { \Delta t}{2} u_{tt}
+ ... Higher Order Terms ... \f]
which again shows first order convergence.

\remarks This Taylor Series was computed using sage, http::www.sagemath.org.  Sage is a Python math program
with symbolic calculus capabilities.  The code follows
\code{.py}
#Euler Explicit Burgers' Equation Finite Difference Form
var('x,t,Dx,Dt,nu')
U = function('U',x,t)
T = taylor(U,( x, 0),(t,0), 3)
S = taylor(U^2,(x, 0),(t,0), 3)
rem = (T.substitute(x=0,t=0) -T.substitute(x=0,t=Dt))/Dt - ( S.substitute(x=0,t=0) - S.substitute(x=-Dx,t=0))/Dx/2
rem += nu*(T.substitute(x=-Dx,t=0) - 2*T.substitute(x=0,t=0) + T.substitute(x=Dx,t=0))/Dx/Dx
print rem.expand()
\endcode

The finite difference nulimit option is \f$ \nu_{reduced} = max( 0.0, \nu - \frac{|u_i| \Delta x}{2} ) \f$.

The finite volume form is implimented with upwinding from the flux through the faces (or mid points)
based on the sign of \f$ \overline{u_{i+1/2}} = \frac{ u_i + u_{i+1} }{2} \f$.

\f$ F_{i+1/2} = \overline{u_{i+1/2}} \frac{u_i}{2} + \nu \frac{ u_{i} - u_{i+1} }{\Delta x}  \f$ - 
if \f$ \overline{u_{i+1/2}} > 0 \f$ ,

\f$ F_{i+1/2} = \overline{u_{i+1/2}} \frac{u_{i+1}}{2} + \nu \frac{ u_{i} - u_{i+1} }{\Delta x}  \f$ -  
if \f$ \overline{u_{i+1/2}} < 0 \f$ and

\f$ F_{i+1/2} =  \nu \frac{ u_{i} - u_{i+1} }{\Delta x}  \f$ - if \f$ \overline{u_{i+1/2}} = 0 \f$ .

Then
\f[ u_i^{n+1} = u_i^n + \frac { \Delta t}{\Delta x} ( F_{i-1/2} - F_{i+1/2} ) \f]
is the full solution equation.

Applying the Taylor Series expansion to all terms for the case \f$ \overline{u_{i+1/2}} > 0 \f$
and \f$ \overline{u_{i-1/2}} > 0 \f$
\f[ u_t + u u_x - \nu u_{xx} = \frac{ \Delta t}{2} u_{tt} + \frac{ \Delta x}{4} ( (u_x)^2 + u u_{xx} ) 
+ ... Higher Order Terms ...\f]
This indicates a slightly lower error for the finite volume form but still only first order convergence.
The sage code for computing the Taylor Series expansions follows:
\code
#Euler Explicit Burgers' Equation Finite Volume Form
var('x,t,Dx,Dt,nu')
U = function('U',x,t)
T = taylor(U,( x, 0),(t,0), 3)
ubm =( T.substitute(x=-Dx,t=0) + T.substitute(x=0,t=0))/2
ubp =( T.substitute(x=Dx,t=0) + T.substitute(x=0,t=0))/2
fp = ubp*T.substitute(x=0,t=0)/2 + nu*(T.substitute(x=0,t=0) - T.substitute(x=Dx,t=0))/Dx
fm = ubm*T.substitute(x=-Dx,t=0)/2 + nu*(T.substitute(x=-Dx,t=0) - T.substitute(x=0,t=0))/Dx
rem =  (T.substitute(x=0,t=0) -T.substitute(x=0,t=Dt))/Dt + ( fm - fp )/Dx
print rem.expand()
\endcode

The finite volume nulimit option is \f$ \nu_{reduced} = max( 0.0, \nu - \frac{|\overline{u_{i+1/2}}| \Delta x}{2} ) \f$.


Summary
-------

The Euler Explicit Method is the simplest method for numerical solution of 
ordinary or partial differential equations. However, it is not very accurate and has
limitted stability. This method is good for demonstrating basic properties of numerical
methods for partial differential equations, such as convergence and stability.
 
*/

class EEWidget : public SolvWidget {

public:
	EEWidget();
	EEWidget ( const EEWidget& other );
	virtual ~EEWidget();
	virtual EEWidget& operator= ( const EEWidget& other );
	virtual bool operator== ( const EEWidget& other ) const;
	
	/** Compute nStep time steps using the current parameters
	@param nStep the number of steps to computer
	*/
	virtual void step ( const size_t nStep );
	
	/// todo - remove from hear - the default should work after fixing step to use the new variables
	virtual void setSize ( const size_t size = 100 );
	
	/** indicate whether equation number equ can be solved in this class
	*/
	virtual bool canSolve(int equ);

protected:
  double *f;
  
  double upwind; /// degree of upwinding for convection, 0 = central, 1.0 = 1st order upwind
  bool nulimit; /// if true, reduce nu to account for upwind effective viscosity
  bool finVol; /// if true, use finite volume form of Burgers' Equation, \f$ u \, u_x \f$ vs \f$ 1/2 (u^2)_x \f$
};

#endif // EEWIDGET_H
