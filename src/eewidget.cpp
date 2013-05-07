/*
    Euler Explicit Method
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

#include "eewidget.h"

void EEWidget::step ( const size_t nStep )
{
    if(N_==0) return;
    if(unstable) return;
    samset = false;
    switch(equation) {
    default:
    case 0:
    case 1:
        advection( nStep);
        return;
    case 2:
    case 3:
        burger( nStep);
        return;
    }
}

void EEWidget::setSize ( const size_t size )
{
    cStep = 0;
    if( size == N_ ) return;
    resize(size);
}

EEWidget::EEWidget(): SolvWidget()
{
    setTitle(tr("Explict Euler 1"));
    plotNameEdit->setText(title);
    setColor(Qt::red);
    unstable = false;
    fluxNames.append(tr("Flux"));
    upwind = 1.0;
    speed_ = 1.0;
    nulimit = false;
    finVol = true;
}

EEWidget::EEWidget ( const EEWidget& other )
{
}

EEWidget::~EEWidget()
{

}

EEWidget& EEWidget::operator= ( const EEWidget& other )
{
    return *this;
}

bool EEWidget::operator== ( const EEWidget& other ) const
{
    return (this == &other);
}

bool EEWidget::canSolve(int equ) {
    return (equ == 0);
}

void EEWidget::advection(const int nStep) {
    double *x = data.value(tr("X"));
    double *u = data.value(tr("U"));
    double *du = data.value(tr("dU"));
    double *e;
    double *d;
    double *f;
    double absc;
    double vred;
    double dx;
    double speed;
    int np;
    int nm;
    if( finVol ) {
        f=fluxes.value(tr("Flux"));
        for( int i=0; i<nStep; i++ ) {
            np = 1;
            for(size_t n=0; n<N_; n++ ) {
                dx = x[np]-x[n];
                if( dx < 0.0) { // either nm<0 or np>=N
                    dx += (xmax_ - xmin_);
                }
                absc = upwind*fabs(speed_)*dx/2.0;
                if(nulimit) {
                    vred = (visc_ > absc )? visc_ : absc;
                } else {
                    vred = visc_ + absc;
                }
                f[n] = speed_*(u[n] + u[np])/2.0 + vred*(u[n] - u[np])/dx ;
                np++;
                if(np == N_) np=0;
            }
            nm = N_-1;
            np = 1;
            for(size_t n=0; n<N_; n++ ) {
                dx = (x[np]-x[nm])/2.0;
                if( dx < 0.0) { // either nm<0 or np>N
                    dx += (xmax_ - xmin_)/2.;
                }
                du[n] = dt/dx*(f[nm] - f[n]);
                nm++;
                np++;
                if( nm == N_ ) nm = 0;
                if( np == N_ ) np = 0;
            }
            for(size_t n=0; n<N_; n++ ) {
                u[n] += du[n];
            }
        }
    } else {
        e = data.value(tr("E"));
        d = data.value(tr("D"));
        for( int i=0; i<nStep; i++ ) {
            nm = N_-1;
            np = 1;
            for(size_t n=0; n<N_; n++ ) {

                dx = (x[np]-x[nm])/2.0;
                if( dx < 0.0) { // either nm<0 or np<N
                    dx += (xmax_ - xmin_)/2.;
                }

                absc = upwind*fabs(speed_)*dx/2.0;
                if(nulimit) {
                    vred = (visc_ > absc )? visc_ : absc;
                } else {
                    vred = visc_ + absc;
                }
                e[n] = speed_*(u[np]-u[nm])/dx/2.0;
                d[n] = vred*(u[nm]-2*u[n]+u[np])/dx/dx;
                du[n] = -e[n] + d[n];
                nm++;
                np++;
                if(nm == N_) nm=0;
                if(np == N_) np=0;
            }
            for(size_t n=0; n<N_; n++ ) {
                u[n] += dt*du[n];
            }
        }
    }
}

void EEWidget::burger(const int nStep) {
    double *x = data.value(tr("X"));
    double *u = data.value(tr("U"));
    double *e;
    double *d;
    double *f;
    double absc;
    double vred;
    double dx;
    double speed;
    size_t np;
    size_t nm;
    if( finVol ) {
        f=fluxes.value(tr("Flux"));
        for( int i=0; i<nStep; i++ ) {
            dx = x[0]-x[N_-1];
            speed = (u[0] + u[N_-1])/2.0;
            absc = upwind*fabs(speed)*dx/2.0;
            if(nulimit) {
                vred = (visc_ > absc )? visc_ : absc;
            } else {
                vred = visc_ + absc;
            }
            vred /= dx;
            for(size_t n=0; n<N_; n++ ) {
                np = n+1;
                if(np == N_) np=0;
                dx = x[np]-x[n];
                speed = (u[np] + u[n])/2.0;
                absc = upwind*fabs(speed)*dx/2.0;
                vred = visc_;
                if(nulimit) {
                    vred = (visc_ > absc )? visc_ : absc;
                } else {
                    vred = visc_ + absc;
                }
                f[n] = speed*(u[n] + u[np])/2.0 + vred*(u[n] - u[np]) ;
            }
            nm = N_-1;
            np = 1;
            for(size_t n=0; n<N_; n++ ) {
                dx = (x[np]-x[nm])/2.0;
                u[n] += dt/dx*(f[n] - f[nm]);
                nm++;
                np++;
                if( nm == N_ ) nm = 0;
                if( np == N_ ) np = 0;
            }
        }
    } else {
        e = data.value(tr("E"));
        d = data.value(tr("D"));
        for( int i=0; i<nStep; i++ ) {
            nm = N_-1;
            np = 1;
            for(size_t n=0; n<N_; n++ ) {
                dx = (x[np]-x[nm])/2.0;
                absc = upwind*fabs(u[n])*dx/2.0;
                if(nulimit) {
                    vred = (visc_ > absc )? visc_ : absc;
                } else {
                    vred = visc_ + absc;
                }
                e[n] = (u[np]*u[np]-u[nm]*u[nm])/4.0;
                d[n] = vred*(u[nm]-2*u[n]+u[np])/dx/dx;
                u[n] += -e[n] + d[n];
                nm++;
                np++;
                if(nm = N_) nm=0;
                if(np = N_) np=0;
            }
        }
    }
}

