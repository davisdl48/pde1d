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


#ifndef IMPWIDGET_H
#define IMPWIDGET_H

#include "solvwidget.h"
#include "myinputs.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <slu_ddefs.h>
#include <supermatrix.h>
//#include <superlu/slu_ddefs.h>
//#include "../superlu/SRC/slu_ddefs.h"

class ImpWidget : public SolvWidget
{

    Q_OBJECT

public:
    ImpWidget ( QWidget* parent = 0 );
    ImpWidget ( const ImpWidget& other );
    virtual ~ImpWidget();
    virtual ImpWidget& operator= ( const ImpWidget& other );
    virtual bool operator== ( const ImpWidget& other ) const;

    virtual void step ( const size_t nStep );
    virtual void setSize ( const size_t size = 100 );
    virtual void setCFL ( const double value = 1.0 ) ;
    const double getImplicit();
    virtual void initSin ( const double value );
    virtual double* getU();
    bool canSolve(int equ) ;
    
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
    double *coef; // coefficient matrix index ( variable + nvar*ntime + nvar*ntime*block )
    int ncoef; // number of coefficients = nblock*nvar*ntime -- used to track changes
    int nblock; //block size
    int ivar; // index of first diagonal element
    int nvar; // number of bases interacting with the diagonal element
    int ntime; // 2 or 3
    int nup; // number of upper blocks or coefficients
    int ndn; // number of lower coefficients or blocks

    double impl;
    double beta;
    double back;
    int ibase;
    double *Ub;
    bool aexist; // check is A needs to be deleted
    bool dirty; // value to check for changes
    char           equed[1];
    yes_no_t       equil;
    trans_t        trans;
    SuperMatrix    A, L, Up;
    SuperMatrix    B, X;
    double         *a;
    int            *asub, *xa;
    int            *perm_c; /* column permutation vector */
    int            *perm_r; /* row permutations from partial pivoting */
    int            *etree;
    void           *work;
    int            info, lwork, nrhs;
    int            i, m, n, nnz;
    double         *rhsb, *rhsx;
    double         *R, *C;
    double         *ferr, *berr;
    double         u, rpg, rcond;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    
    bool transform;
    double *Utran;
    gsl_vector * TranVec; // transform values
    gsl_vector * UVec;  // Physical Space
    gsl_matrix * Mforw; // forward transform matrix = LU of Mback
    gsl_matrix * Mback; // reverse transform matrix
    int matSize_gsl; // flag for control of gsl_matrix allocation
    gsl_permutation * permut;
    int signum;
    
    void updateCoef( int value );
    void fillA();
    void fillB();
    void blockFillA();
    void lspc();
    void lspa();
    void femc();
    void fema() ;
    void reversTrans(double *from,double *to);
    void forwardTrans(double *from,double *to);
    void setupTrans() ;
    
};

#endif // LEASTSQRWIDGET2_H
