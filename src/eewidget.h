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

The derivative or slope of a function <em>u</em> at a location,\f$ x_i \f$ , and time,
\f$ t_n \f$ , with respect to time is defined as 
\f[ \frac{ \partial u }{ \partial t } = \lim_{\delta t \to 0 } 
\frac { u(x_i,t_n + \delta t ) - u(x_i , t_n ) }{\delta t} \f]
where \f$ x_i \f$ is the i'th x point and \f$ t_n \f$ is the n'th time step.  By replacing 
\f$ \delta t \f$ with by a small \f$ \Delta t \f$ on the right side and solving for 
\f$ u(x_i,t_n + \Delta t ) \f$ an approximation for the value of <em>u</em> at the new time
can be obtained.
\f[ u(x_i,t_n + \Delta t ) = u(x_i,t_{n}) + \Delta t \frac{ \partial u }{ \partial t } \f]

\par Time Component 

The Euler Explicit method uses the partial derivative with respect to time at the current
time to predict the value at the next time step, that is
\f[ u_i^{n+1} = u_i^{n} + \Delta t u_t \f]
where \f$ u_i^{n+1} \f$ is an notation for \f$ u(x_i,t_n + \Delta t) \f$ and \f$ u_t \f$
is shorthand notation for \f$ \frac{ \partial u }{ \partial t } \f$.  
 
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
where \f$ |c| \f$ is the absolute value of <em> c </em>.
 
For second derivative terms, a central difference is used, 
 \f$ \nu u_{xx} \approx \nu \frac { u_{i+1} - 2 u_i + u_{i-1} }{ \Delta x^2 } \f$.
   
\par Advection/Diffusion Equation
 
Substituting the above approximations into the Advection/Diffusion equation, 
\f$ u_t + c u_x = \nu u_xx \f$ , results in
\f[ u_i^{n+1} = u_i^{n} + \Delta t \left( -\frac{c}{2 \Delta x} ( u_{i+1} - u_{i-1} ) 
 - \frac{|c|}{2 \Delta x}  ( -u_{i+1} + 2 u_i -u_{i-1} ) 
 + \nu \frac { u_{i+1} - 2 u_i + u_{i-1} }{ \Delta x^2 }
   \right) \f]
or combining a few terms give the difference equation
  \f[ u_i^{n+1} = u_i^{n} + \frac{\Delta t}{\Delta x} 
  \left( -\frac{c}{2} ( u_{i+1} - u_{i-1} ) 
 + \left( \frac{\nu}{\Delta x} + \frac{|c|}{2} \right)  ( u_{i+1} - 2 u_i + u_{i-1} )
   \right) \f].  
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
'volume' average of <em>u</em> , \f$ \bar u \f$ , at the new time level becomes
\f{eqnarray*}{ 
   \overline{u^{n+1}} &=& \frac{1}{\Delta x} \int_{x_{i-1/2}}^{x_{i+1/2}} u(x,t_n + \Delta t)\, dx \\
 &=& \frac {1}{\Delta x} \int_{x_{i-1/2}}^{x_{i+1/2}} u(x,t_n) - \int_{t_n}^{t_n + \Delta t}f_x \,dt \, dx \\
 &=& \frac {1}{\Delta x} \int_{x_{i-1/2}}^{x_{i+1/2}} u(x,t_n)\,dx 
 - \int_{t_n}^{t_n + \Delta t}\int_{x_{i-1/2}}^{x_{i+1/2}}f_x \,dx \, dt
\f}
The integrals can be reversed because we assume that <em>f<sub>x</sub></em> is continuous
and differentiable.  Also, from the Fundamental theorem of calculus (the divergence theorem )
 \f[ \int_{x_{i-1/2}}^{x_{i+1/2}} f_x \, dx = f(x_{i+1/2},t) -f(x_{i-1/2},t) \f] which can be substituted into 
the above equation to yield
 \f[ \overline{u_i^{n+1}} = \overline{ u_i^n} - \frac {1}{\Delta x}
 \int_{t_n}^{t_n + \Delta t} f(x_{i+1/2},t) -f(x_{i-1/2},t) dt \f]
 
Then assuming <em>f</em> is constant with respect to time over the interval from <em>t</em> to 
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
  \f[ u_x|_{i+1/2} = \frac{u_{i+1} - u_i }{\Delta x} \f].  Dropping the overbar for the average of <em>u</em>
gives
  \f[ u_i^{n+1} =  u_i^n + \frac{\Delta t }{\Delta x} \left( \frac{c}{2} ( u_{i-1} - u_{i+1} ) \right)   + \frac{\Delta t }{\Delta x^2} 
  \left( \nu + \frac{|c| \Delta x}{2} \right) (u_{i+1} - 2 u_i + u_{i-1}) \f].
This is the same result as the finite difference equation above.  
 
\par Accuracy
 
To verify the accuracy of this equation substitute the Taylor Series expansion of \f$ u \f$
about \f$ x_i \f$.
 \f[ u_{i-1} = u(x_i - \Delta x,t_n) = u(x_i) - \Delta x u_x(x_i,t_n)
 +  \frac {\Delta x^2}{2} u_{xx}(x_i,t_n) - \frac {\Delta x^3}{3!} u_{xxx}(x_i,t_n) + R_k(x) \f]
 \f[ u_{i+1} = u(x_i + \Delta x,t_n) = u(x_i) + \Delta x u_x(x_i,t_n)
 +  \frac {\Delta x^2}{2} u_{xx}(x_i,t_n) + \frac {\Delta x^3}{3!} u_{xxx}(x_i,t_n) + R_k(x) \f]
where \f$ R_k(X)\f$ is a remander or error term and for functions \f$ u \f$ that are
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
\f[  u_i^n + \Delta t u_t + \frac {\Delta t^2}{2} u_{tt} + \frac {\Delta t^3}{3!} u_{ttt} 
 = u_i^{n} + \frac{\Delta t}{\Delta x}  
    \left(    -\frac{c}{2} ( 2 \Delta x u_x + 2 \frac { \Delta x^3}{3!} u_{xxx} ) 
 + \left( \frac{\nu}{\Delta x} + \frac{|c|}{2} \right)  
 ( \Delta x^2 u_{xx} + 2 \frac {\Delta x^4}{4!} u_{x^{(4)}} )
    \right) + ... \f] .
Collecting some terms, keeping only terms up to first order deltas gives
    \f[  u_t + c u_x - \nu u_{xx}  = \frac{|c| \Delta x}{2} u_{xx} - \frac {\Delta t}{2} u_{tt} \f]
which approaches the advection/diffusion equation as \f$ \Delta x \f$ and \f$ \Delta t \f$ approach zero.
Because the error term on the right hand side involves first order deltas the method is first order.
 
The 'nulimit' option reduces the viscosity by the aparent viscosity do to upwinding through the absolute value of wave speed, <em>|c|</em>.
 \f[ \nu_{reduced} = max( 0.0, \nu - \frac{|c| \Delta x}{2} ) \f]
    
\par Stability
 
For stability analysis consider whether a small error will decay or increase.  
 
\par Burgers' Equation
 
The Euler Explicit Method is the simplest method for numerical solution of 
ordinary or partial differential equations, however it is not very accurate and has
limitted stability.
 
*/

class EEWidget : public SolvWidget {

public:
	EEWidget();
	EEWidget ( const EEWidget& other );
	virtual ~EEWidget();
	virtual EEWidget& operator= ( const EEWidget& other );
	virtual bool operator== ( const EEWidget& other ) const;
	
	virtual void step ( const size_t nStep );
	virtual void setSize ( const size_t size = 100 );
	virtual bool canSolve(int equ);

protected:
  double *f;
  bool nulimit;
  bool finVol;
};

#endif // EEWIDGET_H
