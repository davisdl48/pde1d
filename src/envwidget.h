/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2013  Doug Davis <email>

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


#ifndef ENVWIDGET_H
#define ENVWIDGET_H

#include "solvwidget.h"
#include "myinputs.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
//#include "../superlu/SRC/slu_ddefs.h"

class EnvWidget : public SolvWidget
{

    Q_OBJECT

public:
    EnvWidget ( QWidget* parent = 0 );
    EnvWidget ( const EnvWidget& other );
    virtual ~EnvWidget();
    virtual EnvWidget& operator= ( const EnvWidget& other );
    virtual bool operator== ( const EnvWidget& other ) const;

    virtual void step ( const size_t nStep );
    virtual void setSize ( const size_t size = 100 );
    virtual void setCFL ( const double value = 1.0 ) ;
    const double getImplicit();
    virtual void initSin ( const double value );
    virtual double* getU();
    
public slots:
    void setImplicit ( double value = 5/12.0 ) ;
    void setBackward ( double value = -1/12.0 );
    void setBasis ( int index );
    void setMethod ( int index ) ;


protected:
    QLabel *implLabel;
    MyDoubInput *implInput;
    QLabel *backLabel;
    MyDoubInput *backInput;
    QLabel *weightLabel;
    QComboBox *weightBox;
    QLabel *methodLabel;
    QComboBox *methodBox;

    int method;
    double impl;
    double beta;
    double back;
    int ibase;
    bool dirty; // value to check for changes
    
    int ipad; // padding to the sides of the window
    int winwid; // initial window width
    int winoff; // offset step between windows 
    double * weights; // window weights 
    
    size_t nbas; // number of points in the transform and number of basis functions
                 // nbas == winwid + 2*ipad
    int ntime; // 2 - central Time Implicit, 3 - Adams

    double *Ub; // U at n-1 for adams method or Utransform at n-1 
    size_t nub;            // Ub size = N_/winoff*nbas
    double *Usum;
    gsl_vector_view UbView;
    
    gsl_vector * TranVec; // nbas transform values within envelope
    gsl_vector * UVec;  // Physical Space solution of envelope
    gsl_vector * BVec;  // Right Side
    gsl_vector * CVec;  // far right side for Adams method
    
    gsl_matrix * Left; // Implicit part of pde operator - LU factored
    gsl_permutation * lpermut;
    gsl_matrix * Right;  // Right side
    gsl_matrix * FarRight;  // Right side operator on Ub
    
    gsl_matrix * Mforw; // forward transform matrix = LU of Mback
    gsl_matrix * Mback; // reverse transform matrix
    
   
    gsl_permutation * permut;
    int signum;
    
    void updateCoef( int value );
    void lspc();
    void lspa();
    void femc();
    void fema();
    void rk();
    void setupTrans(int size) ;
    void allocate_gsl(int size);
    
};

#endif // LEASTSQRWIDGET2_H
