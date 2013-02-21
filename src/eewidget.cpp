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
  if(N==0) return;
  for(size_t n = 0; n < nStep; n++) {
    if(CFL > 0) {
      for(size_t i=0; i<N; i++) {
	f[i]=U[i];
      }
    }else{
      for(size_t i=0; i<(N-1); i++) {
	f[i]=U[i+1];       
      }
    }  
    for(size_t i=1; i<N; i++) {
      U[i] = U[i]+CFL*(f[i-1]-f[i]);
      //std::cout << x[i] << '\t' << U[i] << std::endl;
    }
    U[0] = U[0]+CFL*(f[N-1]-f[0]);
    //std::cout << std::endl << std::endl;
    cStep++;
    totCFL += CFL;
  }
  /*
  for(size_t i=0; i<N; i++) {
    fout << x[i] << '\t' << U[i] << std::endl;
  }
  fout << std::endl << std::endl;
  */
  //gnuplot_close(h1);
}

void EEWidget::setSize ( const size_t size )
{

  cStep = 0;
  if( size == N ) return;
  if( N != 0) {
    delete[] U;
    delete[] f;
    delete[] x;
    delete[] ideal;
  }
  N = size;
  U = new double[N];
  f = new double[N];
  x = new double[N];
  ideal = new double[N];
  initSin(cycles);
}

EEWidget::EEWidget(): SolvWidget()
{
  setTitle(tr("Explict Euler 1"));
  plotNameEdit->setText(title);
  setColor(Qt::red);
}

EEWidget::EEWidget ( const EEWidget& other )
{
}

EEWidget::~EEWidget()
{
  if(N != 0) {
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

