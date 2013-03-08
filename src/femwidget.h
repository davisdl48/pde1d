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


#ifndef LS2WIDGET_H
#define LS2WIDGET_H
#include "myinputs.h"
#include <gsl/gsl_vector.h>

#include "solvwidget.h"


class FEMWidget : public SolvWidget {
  Q_OBJECT

public:
	FEMWidget(QWidget *parent=0);
	FEMWidget ( const FEMWidget& other );
	virtual ~FEMWidget();
	virtual FEMWidget& operator= ( const FEMWidget& other );
	virtual bool operator== ( const FEMWidget& other ) const;
	
	virtual void step ( const size_t nStep );
	virtual void setSize ( const size_t size = 100 );
        const double getAlpha();
    
    public slots:
      void setAlpha(double value = 0.5);
      void setBasis(int index);
      
	
protected:
  QLabel *alphaLabel;
  MyDoubInput *alphaInput;
  QLabel *weightLabel;
  QComboBox *weightBox;
  
  double a;
  double b;
  bool cosBas;
  gsl_vector * X;
  gsl_vector * DIAG;
  gsl_vector * E;
  gsl_vector * F;
  gsl_vector * B;
};

#endif // LS2WIDGET_H
