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


#ifndef SIMPIMPWIDGET_H
#define SIMPIMPWIDGET_H

#include "solvwidget.h"
#include "myinputs.h"
#include <gsl/gsl_vector.h>


class SimpImpWidget : public SolvWidget
{
  Q_OBJECT

public:
    SimpImpWidget(QWidget *parent=0);
    SimpImpWidget(const SimpImpWidget& other);
    virtual ~SimpImpWidget();
    virtual SimpImpWidget& operator=(const SimpImpWidget& other);
    virtual bool operator==(const SimpImpWidget& other) const;
    double getUpwind() ;
    double getImpl() ;
    virtual void step(size_t nStep=1);
    virtual void setSize(const size_t value);
    
public slots:
  void setImpl(double value = 0.0) ;
  void setUpwind(double value = 0.5) ;
  
protected:
  double impl;
  double upwind;
  QLabel *impLabel;
  MyDoubInput *impInput;
  QLabel *upwLabel;
  MyDoubInput *upwInput;
  
  gsl_vector * X;
  gsl_vector * DIAG;
  gsl_vector * E;
  gsl_vector * F;
  gsl_vector * B;
};

#endif // SIMPIMPWIDGET_H
