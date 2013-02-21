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


#ifndef LEASTSQRWIDGET_H
#define LEASTSQRWIDGET_H

#include "solvwidget.h"
#include <knuminput.h>
#include <gsl/gsl_vector.h>

class LeastSqrWidget : public SolvWidget
{
    Q_OBJECT

public:
    LeastSqrWidget ( QWidget *parent = 0 );
    LeastSqrWidget ( const LeastSqrWidget& other );
    virtual ~LeastSqrWidget();
    virtual LeastSqrWidget& operator= ( const LeastSqrWidget& other );
    virtual bool operator== ( const LeastSqrWidget& other ) const;

    virtual void step ( size_t nStep = 1 );
    virtual void setSize ( const size_t value );
    const double getAlpha();

public slots:
    void setAlpha ( double value = 0.5 );
    void setBasis ( int index );


protected:
    QLabel *alphaLabel;
    KDoubleNumInput *alphaInput;
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

#endif // LEASTSQRWIDGET_H
