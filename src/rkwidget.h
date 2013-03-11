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


#ifndef RKWIDGET_H
#define RKWIDGET_H

#include "solvwidget.h"
#include "myinputs.h"
#include <superlu/slu_ddefs.h>
//#include "../superlu/SRC/slu_ddefs.h"

class RKWidget : public SolvWidget
{

    Q_OBJECT

public:
    RKWidget ( QWidget* parent = 0 );
    RKWidget ( const RKWidget& other );
    virtual ~RKWidget();
    virtual RKWidget& operator= ( const RKWidget& other );
    virtual bool operator== ( const RKWidget& other ) const;

    virtual void step ( const size_t nStep );
    virtual void setSize ( const size_t size = 100 );
    virtual void setCFL ( const double value = 1.0 ) ;
    const double getImplicit();
    bool canSolve(int equ) ;

public slots:
    void setBasis ( int index );
    void setMethod ( int index ) ;


protected:
    QLabel *weightLabel;
    QComboBox *weightBox;
    QLabel *methodLabel;
    QComboBox *methodBox;

    int method;
    double *coef; // coefficient matrix index ( variable + nvar*2 + nvar*2*block )
    int ncoef; // number of coefficients = nblock*nvar*2 -- used to track changes
    int nblock; //block size
    int ivar; // index of first diagonal element
    int nvar; // number of bases interacting with the diagonal element
    int nup; // number of upper blocks or coefficients
    int ndn; // number of lower coefficients or blocks

    // Butcher Tablue Values for RK integrations
    double *b_a;// assume lower triangular
    double *b_b;
    double *b_c;
    double *b_k;// dimension N_*nstage*n_b_k
    int nStage;
    int n_b_k; // for RKF adaptive methods - modified Butcher Tablue size

    int ibase;
    bool aexist; // check if A needs to be deleted
    //bool dirty; // value to indicate changes that require a new A matrix
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

    void updateCoef( int value );
    void fillA();
    void fillB(int stage);
    void blockFillA();
    void setNStage(int arg1);
};

#endif // LEASTSQRWIDGET2_H
