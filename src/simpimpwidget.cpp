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
    impl = 0.5;
    upwind = 0.0;
    impLabel = new QLabel(QString::fromLocal8Bit("Implicit"));
    impInput = new KDoubleNumInput(0.0,1.0,0.5,parent,0.01,6);
   
    connect(impInput,SIGNAL(valueChanged(double)),this,SLOT(setImpl(double)));
    upwLabel = new QLabel(QString::fromLocal8Bit("Upwinding"));
    upwInput = new KDoubleNumInput(-0.1,1.1,0.0,parent,0.01,6);
    connect(upwInput,SIGNAL(valueChanged(double)),this,SLOT(setUpwind(double)));
    verticalLayout->insertWidget(4,impLabel); 
    verticalLayout->insertWidget(5,impInput);
    verticalLayout->insertWidget(6,upwLabel);
    verticalLayout->insertWidget(7,upwInput);
    setColor( Qt::green );

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
    impl = value;

}

void SimpImpWidget::setUpwind(double value) {
    upwind = value;
}

void SimpImpWidget::setSize(const size_t value)
{
  //std::cout << "SimpImpLW::setSize( " << value << " )\n";
  if( value == N) return;
  cStep = 0;
  if( N != 0) {
    gsl_vector_free (DIAG);
    gsl_vector_free (E);// upper 
    gsl_vector_free (F);// lower
    gsl_vector_free (B);
    gsl_vector_free (X);
    delete[] U;
    delete[] x;
    delete[] ideal;
  }
  N = value;
  U = new double[N];
  x = new double[N];
  ideal = new double[N];
  DIAG = gsl_vector_alloc (N);
  E = gsl_vector_alloc (N);
  F = gsl_vector_alloc (N);
  B = gsl_vector_alloc (N);
  X = gsl_vector_alloc (N);
  initSin(cycles);
}

void SimpImpWidget::step(size_t nStep) {
  if(N==0) return;
  double f0,fm,fp;
  double ufun;
  if(CFL > 0 ) {
    fp= (upwind-1)*CFL/2;
    f0 = -upwind*CFL;
    fm = (1+upwind)*CFL/2;
  }else{
    fp = -(1+upwind)*CFL/2;
    f0 = upwind*CFL;
    fm = (1-upwind)*CFL/2;
  }
  for(size_t n = 0; n< nStep; n++) {
    ufun = fm*U[N-1]+f0*U[0]+fp*U[1];
    gsl_vector_set (B,0,ufun);
    for(size_t i=1; i<(N-1); i++) {
      ufun = fm*U[i-1]+f0*U[i]+fp*U[i+1];
      gsl_vector_set (B,i,ufun);
    }
    ufun = fm*U[N-2]+f0*U[N-1]+fp*U[0];
    gsl_vector_set (B,N-1,ufun);
    gsl_vector_set_all(DIAG,1-f0*impl);
    gsl_vector_set_all(E,-fp*impl);
    gsl_vector_set_all(F,-fm*impl);
    gsl_linalg_solve_cyc_tridiag (DIAG,E,F,B,X);
    for(size_t i=0; i<N; i++) {
      U[i] = U[i]+gsl_vector_get(X,i);
      //std::cout << x[i] << '\t' << U[i] << std::endl;
    }
    //std::cout << std::endl << std::endl;
    cStep++;
    totCFL += CFL;
  }
  return;
  /*
  for(size_t i=0; i<N; i++) {
    fout << x[i] << '\t' << U[i] << std::endl;
  }
  fout << std::endl << std::endl;
  */
  //gnuplot_close(h1);
}



