/*
    Euler Explicit Method with posibility to add first order upwinding
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
  for(size_t n = 0; n < nStep; n++) {
    if(CFL > 0) {
      for(size_t i=0; i<N_; i++) {
	f[i]=U_[i];
      }
    }else{
      for(size_t i=0; i<(N_-1); i++) {
	f[i]=U_[i+1];       
      }
    }  
    for(size_t i=1; i<N_; i++) {
      U_[i] = U_[i]+CFL*(f[i-1]-f[i]);
      if(U_[i] > 1e16) unstable = true;
      //std::cout << X_[i] << '\t' << U_[i] << std::endl;
    }
    U_[0] = U_[0]+CFL*(f[N_-1]-f[0]);
    //std::cout << std::endl << std::endl;
    cStep++;
    totCFL += CFL;
  }
  /*
  for(size_t i=0; i<N_; i++) {
    fout << X_[i] << '\t' << U_[i] << std::endl;
  }
  fout << std::endl << std::endl;
  */
  //gnuplot_close(h1);
}

void EEWidget::setSize ( const size_t size )
{

  cStep = 0;
  if( size == N_ ) return;
  if( N_ != 0) {
    delete[] U_;
    delete[] f;
    delete[] X_;
    delete[] Ideal_;
  }
  N_ = size;
  U_ = new double[N_];
  f = new double[N_];
  X_ = new double[N_];
  Ideal_ = new double[N_];
  initSin(cycles);
}

EEWidget::EEWidget(): SolvWidget()
{
  setTitle(tr("Explict Euler 1"));
  plotNameEdit->setText(title);
  setColor(Qt::red);
  unstable = false;
}

EEWidget::EEWidget ( const EEWidget& other )
{
}

EEWidget::~EEWidget()
{
  if(N_ != 0) {
    delete[] f;
  }
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

