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


#include "leastsqrwidget.h"
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <iomanip>

LeastSqrWidget::LeastSqrWidget ( QWidget *parent ) : SolvWidget ( parent )
{
    setTitle(tr("Least Squares"));
    plotNameEdit->setText ( title );
    alphaLabel = new QLabel ( QString::fromLocal8Bit ( "Alpha" ) );
    alphaInput = new MyDoubInput ( 0.5, this, 0.0, 1.0, 0.01, 6 );
    connect ( alphaInput, SIGNAL ( valueChanged ( double ) ), this, SLOT ( setAlpha ( double ) ) );
    verticalLayout->insertWidget ( 4, alphaLabel );
    verticalLayout->insertWidget ( 5, alphaInput );
    weightLabel = new QLabel ( tr ( "Basis Functions" ) );
    weightBox = new QComboBox ( this );
    weightBox->addItem ( tr ( "Linear Basis" ) );
    weightBox->addItem ( tr ( "Cosine Basis" ) );
    weightBox->setCurrentIndex ( 0 );
    verticalLayout->insertWidget ( 6, weightLabel );
    verticalLayout->insertWidget ( 7, weightBox );
    connect ( weightBox, SIGNAL ( activated ( int ) ), this, SLOT ( setBasis ( int ) ) );
    setColor ( Qt::blue );
    a=-1;
    setAlpha (0.5);
    setBasis(0);
  unstable = false;
}

LeastSqrWidget::LeastSqrWidget ( const LeastSqrWidget& other )
{

}

LeastSqrWidget::~LeastSqrWidget()
{
    if ( N_ != 0 ) {
        gsl_vector_free ( DIAG );
        gsl_vector_free ( E ); // upper
        gsl_vector_free ( F ); // lower
        gsl_vector_free ( B );
        gsl_vector_free ( X );
    }
}

LeastSqrWidget& LeastSqrWidget::operator= ( const LeastSqrWidget& other )
{
    return *this;
}

bool LeastSqrWidget::operator== ( const LeastSqrWidget& other ) const
{
    return (this == &other);
}

void LeastSqrWidget::setAlpha ( double value )
{ 
    b = 1 - value;
    if( value == a) return;
    a = value;
    alphaInput->setValue(a);
}

const double LeastSqrWidget::getAlpha()
{
    return a;
}

void LeastSqrWidget::setBasis ( int index )
{
    switch ( index ) {
    case 0:
        cosBas = false;
        break;
    case 1:
        cosBas = true;
        break;
    default:
        cosBas = false;
        weightBox->setCurrentIndex ( 0 );
    }
}


void LeastSqrWidget::setSize ( const size_t value )
{
    std::cout << "LeastSqrLW::setSize( " << value << " )\n";
    if ( value == N_ ) return;
    cStep = 0;
    if ( N_ != 0 ) {
        gsl_vector_free ( DIAG );
        gsl_vector_free ( E ); // upper
        gsl_vector_free ( F ); // lower
        gsl_vector_free ( B );
        gsl_vector_free ( X );
        delete[] U_;
        delete[] X_;
        delete[] Ideal_;
    }
    N_ = value;
    U_ = new double[N_];
    X_ = new double[N_];
    Ideal_ = new double[N_];
    DIAG = gsl_vector_alloc ( N_ );
    E = gsl_vector_alloc ( N_ );
    F = gsl_vector_alloc ( N_ );
    B = gsl_vector_alloc ( N_ );
    X = gsl_vector_alloc ( N_ );
    initSin ( cycles );
}

void LeastSqrWidget::step ( size_t nStep )
{
    if ( N_ == 0 ) return;
    double f00, f01, um1, u00, up1;
    double b2, ab, p2;
    double ufun;
    if(unstable) return;
    //gnuplot_ctrl * h1;
    //h1 = gnuplot_init();
    b2 = b * b * CFL * CFL;
    ab = a * b * CFL * CFL;
    //cos basis test
    if(cosBas) {
        p2 = pi * pi;
        f00 = ( 3 + b2 * p2 ) / 4.;
        f01 = ( 1.0 - b2 * p2 ) / 8.;
        um1 = ( 1 + 4 * CFL + ab * p2 ) / 8.0;
        u00 = ( 3.0 - ab * p2 ) / 4.;
        up1 = ( 1 - 4 * CFL + ab * p2 ) / 8.0;
    } else {
        // linear basis
        f00 = 2.0 * ( 1.0 / 3.0 + b2 );
        f01 = 1.0 / 6.0 - b2;
        um1 = ( 1.0 / 3.0 + CFL ) / 2.0 + ab;
        u00 = 2 * ( 1.0 / 3.0 - ab );
        up1 = ( 1.0 / 3.0 - CFL ) / 2.0 + ab;
    }
    for ( size_t n = 0; n < nStep; n++ ) {
        ufun = um1 * U_[N_ - 1] + u00 * U_[0] + up1 * U_[1];
        gsl_vector_set ( B, 0, ufun );
        for ( size_t i = 1; i < ( N_ - 1 ); i++ ) {
            ufun = um1 * U_[i - 1] + u00 * U_[i] + up1 * U_[i + 1];
            gsl_vector_set ( B, i, ufun );
        }
        ufun = um1 * U_[N_ - 2] + u00 * U_[N_ - 1] + up1 * U_[0];
        gsl_vector_set ( B, N_ - 1, ufun );
        gsl_vector_set_all ( DIAG, f00 );
        gsl_vector_set_all ( E, f01 );
        gsl_linalg_solve_symm_cyc_tridiag ( DIAG, E, B, X );
        for ( size_t i = 0; i < N_; i++ ) {
            U_[i] = gsl_vector_get ( X, i );
	    if(U_[i] > 1e16) unstable = true;
            //std::cout << X_[i] << '\t' << U_[i] << std::endl;
        }
        //std::cout << std::endl << std::endl;
        cStep++;
        totCFL += CFL;
    }
    return;/*
  for(size_t i=0; i<N_; i++) {
    fout << X_[i] << '\t' << U_[i] << std::endl;
  }
  fout << std::endl << std::endl;
  */
    //gnuplot_close(h1);
}
bool LeastSqrWidget::canSolve(int equ) {
    return (equ == 0);
}
