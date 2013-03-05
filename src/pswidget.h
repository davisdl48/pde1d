/*
    Pseudo Spectral solver for U_t + e_*E_(U)_x  + d_*D_(U)_xx = 0.0
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


#ifndef PSWIDGET_H
#define PSWIDGET_H

#include "solvwidget.h"
#include <knuminput.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_errno.h>
//#include "../superlu/SRC/slu_ddefs.h"

class PSWidget : public SolvWidget
{

    Q_OBJECT

public:
    PSWidget ( QWidget* parent = 0 );
    PSWidget ( const PSWidget& other );
    virtual ~PSWidget();
    virtual PSWidget& operator= ( const PSWidget& other );
    virtual bool operator== ( const PSWidget& other ) const;

    virtual void step ( const size_t nStep );
    virtual void setSize ( const size_t size = 100 );
    virtual void setCFL ( const double value = 1.0 ) ;
    const double getImplicit();
    virtual void initSin ( const double value );
    virtual double* getU();
    void setNStage(int arg1) ;

public slots:
    void setImplicit ( double value = 5/12.0 ) ;
    void setBackward ( double value = -1/12.0 );
    void setMethod ( int index ) ;
	  

protected:
    QLabel *implLabel;
    KDoubleNumInput *implInput;
    QLabel *backLabel;
    KDoubleNumInput *backInput;
    QLabel *weightLabel;
    QComboBox *weightBox;
    QLabel *methodLabel;
    QComboBox *methodBox;
    QLabel *pdeLabel;
    QComboBox *pdeBox;

    int method;
    double impl;
    double beta;
    double back;

    int ntime; // 2 - central Time Implicit, 3 - Adams

    gsl_fft_real_wavetable * real_g;
    gsl_fft_halfcomplex_wavetable * hc_g;
    gsl_fft_real_workspace * work_g;

    double **data; // array of pointers to double
    size_t narrays; // number of data arrays of size N_
    /* data[0] -> time n+1 data
     * data[1] -> time n data
     * data[2] -> time n-1 data
     * data[3 to 3+nk] -> Runge Kutta stage data
     */
    size_t nrk; // number of Runge Kutta Stages
    size_t dsize; // allocated size of data arrays = N_
    
    

    void allocateData(size_t rksize) ;
    void freeData();

    void lspc();
    void lspa();
    void femc();
    void fema();
    void rk();
    void setupTrans();
    void phaser() ;

    // Butcher Tablue Values for RK integrations
    double *b_a;// assume lower triangular
    double *b_b;
    double *b_c;
    int nStage;
    int n_b_k; // for RKF adaptive methods - modified Butcher Tablue size


};

#endif // LEASTSQRWIDGET2_H
