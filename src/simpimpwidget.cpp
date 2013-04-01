/*
    <one line to give the library's name and an idea of what it does.>
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


#include "simpimpwidget.h"
#include <iostream>
#include <gsl/gsl_linalg.h>


SimpImpWidget::SimpImpWidget(QWidget *parent): SolvWidget(parent)
{
    setTitle(tr("Simple Implicit"));
    plotNameEdit->setText(title);
    impLabel = new QLabel(QString::fromLocal8Bit("Implicit"));
    impInput = new MyDoubInput(0.5,this,0.0,1.0,0.01,6);

    connect(impInput,SIGNAL(valueChanged(double)),this,SLOT(setImpl(double)));
    upwLabel = new QLabel(QString::fromLocal8Bit("Upwinding"));
    upwInput = new MyDoubInput(0.0,this,-0.1,1.1,0.01,6);
    connect(upwInput,SIGNAL(valueChanged(double)),this,SLOT(setUpwind(double)));
    verticalLayout->insertWidget(4,impLabel);
    verticalLayout->insertWidget(5,impInput);
    verticalLayout->insertWidget(6,upwLabel);
    verticalLayout->insertWidget(7,upwInput);
    setColor( Qt::green );
    unstable = false;
    impl = -1; // to prevent auto return
    setImpl(0.5);
    upwind = -1; // to prevent auto return
    setUpwind(0.0);

}

SimpImpWidget::SimpImpWidget(const SimpImpWidget& other)
{

}

SimpImpWidget::~SimpImpWidget()
{
    gsl_vector_free (DIAG);
    gsl_vector_free (E);// upper
    gsl_vector_free (F);// lower
    gsl_vector_free (B);
    gsl_vector_free (X);

}

SimpImpWidget& SimpImpWidget::operator=(const SimpImpWidget& other)
{
    return *this;
}

bool SimpImpWidget::operator==(const SimpImpWidget& other) const
{
    return (this == &other);
}

void SimpImpWidget::setImpl(double value) {
    if(impl == value) return;
    impl = value;
    impInput->setValue(impl);

}

void SimpImpWidget::setUpwind(double value) {
  if(upwind == value) return;
  upwind = value;
  upwInput->setValue(upwind);
}

void SimpImpWidget::setSize(const size_t value)
{
    //std::cout << "SimpImpLW::setSize( " << value << " )\n";
    if( value == N_) return;
    if( N_ != 0) { // free derived class storage
        gsl_vector_free (DIAG);
        gsl_vector_free (E);// upper
        gsl_vector_free (F);// lower
        gsl_vector_free (B);
        gsl_vector_free (X);
    }

    // allocate local storage
    DIAG = gsl_vector_alloc (value);
    E = gsl_vector_alloc (value);
    F = gsl_vector_alloc (value);
    B = gsl_vector_alloc (value);
    X = gsl_vector_alloc (value);
    resize(value); // free and reallocate SolvWidget storage and initialize
}

void SimpImpWidget::step(size_t nStep) {
    if(N_==0) return;
    double f0,fm,fp;
    double v0,vm,vp;
    int nm,np;
    double ufun;
    if(unstable) return;
    for(size_t n = 0; n< nStep; n++) {
        Efunc(U_);
        Dfunc(U_);
        for(size_t nn=0; nn<N_; nn++) {
            nm = nn-1;
            nm = (nm < 0 )? N_+nm : nm;
            np = nn+1;
            np = (np >= N_)? np-N_ : np;
            if(CFL*J_[nn] > 0 ) {
                fp= e_*(upwind-1)*CFL/2;
                f0 = -e_*upwind*CFL;
                fm = e_*(1+upwind)*CFL/2;
            } else {
                fp = -e_*(1+upwind)*CFL/2;
                f0 = e_*upwind*CFL;
                fm = e_*(1-upwind)*CFL/2;
            }
            vp = d_*visc_*CFL/dx;
            v0 = -2*d_*visc_*CFL/dx;
            vm = d_*visc_*CFL/dx;
            ufun = fm*E_[nm]+f0*E_[nn]+fp*E_[np] + vm*D_[nm]+v0*D_[nn]+vp*D_[np];
            //std::cout <<  ufun << '\t' << U_[nn] << std::endl;
            gsl_vector_set (B,nn,ufun);
            gsl_vector_set(DIAG,nn,1-(f0*J_[nn]+v0)*impl);
            gsl_vector_set(E,nn,-(fp*J_[np]+vp)*impl);
            gsl_vector_set(F,nn,-(fm*J_[nm]+vm)*impl);
        }
        //std::cout << std::endl << std::endl;
        gsl_linalg_solve_cyc_tridiag (DIAG,E,F,B,X);
        for(size_t i=0; i<N_; i++) {
            //std::cout <<  gsl_vector_get(X,i) << '\t' << U_[i] << std::endl;
            U_[i] = U_[i]+ gsl_vector_get(X,i);
	    if(U_[i] > 1e16) unstable = true;
        }
        //std::cout << std::endl << std::endl;
        cStep++;
        totCFL += CFL;
    }
    return;
    /*
    for(size_t i=0; i<N_; i++) {
      fout << X_[i] << '\t' << U_[i] << std::endl;
    }
    fout << std::endl << std::endl;
    */
    //gnuplot_close(h1);
}
bool SimpImpWidget::canSolve(int equ) {
    return true;
}



