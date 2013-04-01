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


#include "femwidget.h"
#include <gsl/gsl_linalg.h>
#include <iostream>

FEMWidget::FEMWidget(QWidget *parent) : SolvWidget(parent)
{
    setTitle(tr("Finite Element Method"));
    plotNameEdit->setText(title);
    setColor( Qt::magenta );
    alphaLabel = new QLabel(QString::fromLocal8Bit("Implicit"));
    alphaInput = new MyDoubInput(a,this,0.0,1.0,0.01,6);
    verticalLayout->insertWidget(4,alphaLabel);
    verticalLayout->insertWidget(5,alphaInput);
    connect(alphaInput,SIGNAL(valueChanged(double)),this,SLOT(setImplicit(double)));
    weightLabel = new QLabel( tr("Basis Functions"));
    weightBox = new QComboBox(this);
    weightBox->addItem(tr("Linear Basis"));
    weightBox->addItem(tr("Cosine Basis"));
    weightBox->setCurrentIndex ( 0 );
    verticalLayout->insertWidget(6,weightLabel);
    verticalLayout->insertWidget(7,weightBox);
    connect(weightBox,SIGNAL(activated(int)),this,SLOT(setBasis(int)));
    a=-1;
    setImplicit(0.5);
    setBasis(0);
  unstable = false;
}

FEMWidget::FEMWidget ( const FEMWidget& other )
{

}

FEMWidget::~FEMWidget()
{

    if(N_ != 0) {
        gsl_vector_free (DIAG);
        gsl_vector_free (E);// upper
        gsl_vector_free (F);// lower
        gsl_vector_free (B);
        gsl_vector_free (X);
    }
}

FEMWidget& FEMWidget::operator= ( const FEMWidget& other )
{
    return *this;
}

bool FEMWidget::operator== ( const FEMWidget& other ) const
{
    return (this == &other);
}

void FEMWidget::setImplicit(double value) {
  
    b = 1 - value;
    if( value == a) return;
    a = value;
    alphaInput->setValue(a);
}

const double FEMWidget::getImplicit()
{
    return a;
}

void FEMWidget::setBasis ( int index )
{
    switch ( index ) {
    case 0:
        cosBas=false;
        break;
    case 1:
        cosBas=true;
        break;
    default:
        cosBas=false;
        weightBox->setCurrentIndex ( 0 );
    }
}

void FEMWidget::step ( const size_t nStep )
{
    if(N_==0) return;
    double f0,fm,fp,um,u0,up;
    double ufun;
    double tt,six;
    double tf,eith;
    if(unstable) return;
    b = 1-a;
    if(cosBas) {
        //Cosine Bases
        eith= 1.0/8.0;
	tf = 3.0/4.0;
        fm = eith + b*CFL/2;
        f0 = tf;
        fp = eith - b*CFL/2;
        um = eith - a*CFL/2;
        u0 = tf;
        up = eith + a*CFL/2 ;
    } else {
        six=1.0/6.0;
        tt=2.0/3.0;
        fm = six + b*CFL/2;
        f0 = tt;
        fp = six - b*CFL/2;
        um = six - a*CFL/2;
        u0 = tt;
        up = six + a*CFL/2 ;
    }
    for(size_t n = 0; n< nStep; n++) {
        ufun = fm*U_[N_-1]+f0*U_[0]+fp*U_[1];
        gsl_vector_set (B,0,ufun);
        for(size_t i=1; i<(N_-1); i++) {
            ufun = fm*U_[i-1]+f0*U_[i]+fp*U_[i+1];
            gsl_vector_set (B,i,ufun);
        }
        ufun = fm*U_[N_-2]+f0*U_[N_-1]+fp*U_[0];
        gsl_vector_set (B,N_-1,ufun);
        gsl_vector_set_all(DIAG,u0);
        gsl_vector_set_all(E,up);
        gsl_vector_set_all(F,um);
        gsl_linalg_solve_cyc_tridiag (DIAG,E,F,B,X);
        for(size_t i=0; i<N_; i++) {
            U_[i] = gsl_vector_get(X,i);
	    if(U_[i] > 1e16) unstable = true;
            //std::cout << X_[i] << '\t' << U_[i] << std::endl;
        }
        //std::cout << std::endl << std::endl;
        cStep++;
        totCFL += CFL;
    }
    return;

}

void FEMWidget::setSize ( const size_t size )
{
    if( size == N_) return;
    cStep = 0;
    if( N_ != 0) {
        gsl_vector_free (DIAG);
        gsl_vector_free (E);// upper
        gsl_vector_free (F);// lower
        gsl_vector_free (B);
        gsl_vector_free (X);
        delete[] U_;
        delete[] X_;
        delete[] Ideal_;
    }
    N_ = size;
    U_ = new double[N_];
    X_ = new double[N_];
    Ideal_ = new double[N_];
    DIAG = gsl_vector_alloc (N_);
    E = gsl_vector_alloc (N_);
    F = gsl_vector_alloc (N_);
    B = gsl_vector_alloc (N_);
    X = gsl_vector_alloc (N_);
    initSin(cycles);
}
bool FEMWidget::canSolve(int equ) {
    return (equ == 0);
}

