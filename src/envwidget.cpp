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

#include "envwidget.h"
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <cmath>

const bool pdeSolvers[4][5] =  {
                   {true, false, true, false, true},
                   {false, false, false, false, false},
                   {false, false, false, false, false},
                   {false, false, false, false, false} };
		   
EnvWidget::EnvWidget ( QWidget* parent ) : SolvWidget(parent)
{   int iwid = 4;
    setTitle(tr("Envelope Method"));

    methodLabel = new QLabel ( tr ( "Method" ) );
    methodBox = new QComboBox ( this );   
    methodBox->setModel(methodModel);
    methodModel->appendRow( new QStandardItem ( tr("LSP - Cent. time") )); // (0)
    methodModel->appendRow( new QStandardItem ( tr("LSP - Adams") ));      // (1)
    methodModel->appendRow( new QStandardItem ( tr("FEM - Cent. time") )); // (2)
    methodModel->appendRow( new QStandardItem ( tr("FEM - Adams") ));      // (3)
    methodModel->appendRow( new QStandardItem ( tr("Runge Kutta") ));      // (4)
    connect ( methodBox, SIGNAL( activated(int) ), this, SLOT( setMethod(int) ) );
    verticalLayout->insertWidget ( iwid++, methodLabel );
    verticalLayout->insertWidget ( iwid++, methodBox );
    weightLabel = new QLabel ( tr ( "Basis Functions" ) );
    weightBox = new QComboBox ( this );
    weightBox->addItem ( tr ( "Linear by 4" ) );
    weightBox->addItem ( tr ( "Cosine by 8" ) );
    weightBox->addItem ( tr ( "Cosine by 16" ) );
    weightBox->addItem ( tr ( "Cosine by 15" ) );
    verticalLayout->insertWidget ( iwid++, weightLabel );
    verticalLayout->insertWidget ( iwid++, weightBox );
    connect ( weightBox, SIGNAL( activated(int) ), this, SLOT( setBasis(int) ) );
    plotNameEdit->setText ( title );
    implLabel = new QLabel ( tr ( "Implicit" ) );
    implInput = new MyDoubInput ( 5/12.0, this, 0.0, 1.0, 1e-14, 14 );
    implInput->setValue(5/12.0);
    connect ( implInput, SIGNAL ( valueChanged ( double ) ), this, SLOT ( setImplicit ( double ) ) );
    verticalLayout->insertWidget ( iwid++, implLabel );
    verticalLayout->insertWidget ( iwid++, implInput );
    backLabel = new QLabel ( tr ( "Backward" ) );
    backInput = new MyDoubInput ( -1/12.0, this, -1.0, 1.0, 1e-14, 14 );
    backInput->setValue(-1/12.0);
    connect ( backInput, SIGNAL ( valueChanged ( double ) ), this, SLOT ( setBackward(double)) );
    verticalLayout->insertWidget ( iwid++, backLabel );
    verticalLayout->insertWidget ( iwid++, backInput );

    /* Set the default values */
    setColor ( Qt::magenta );
    nbas = 0;

    winwid = 0;
    winoff = 0;


    dirty = true;
    ibase=-1;
    method=-1;
    setBasis(3);
    setMethod(2);
  unstable = false;
    //set_default_options(&options);
}

EnvWidget::EnvWidget ( const EnvWidget& other )
{

}

EnvWidget::~EnvWidget()
{
    if ( N_ != 0 ) {
        if(nbas) {
            delete[] Ub;
        }
    }

    if( nbas ) {
        // deallocate matricies
        gsl_vector_free(UVec);
        gsl_vector_free(BVec);
        gsl_vector_free(CVec);
        gsl_vector_free(TranVec);
        gsl_matrix_free(Mforw);
        gsl_permutation_free(permut);
        gsl_matrix_free(Mback);
        gsl_matrix_free(Left);
        gsl_permutation_free(lpermut);
        gsl_matrix_free(Right);
        gsl_matrix_free(FarRight);
    }
}

EnvWidget& EnvWidget::operator= ( const EnvWidget& other )
{
    return *this;
}

bool EnvWidget::operator== ( const EnvWidget& other ) const
{
    return (this == &other);
}


void EnvWidget::setBasis ( int index )
{
    double dthet;
    if(ibase != index) dirty = true;
    ibase = index;
    switch(ibase) {
    default:
        ibase = 0;
    case 0: // linear by 4
        winwid=4;
        winoff=2;
        weights = new double[winwid];
        weights[0]=0.0;
        weights[1]=weights[3]=0.5;
        weights[2]=1.0;
        setupTrans(winwid);
        break;
    case 1: // raised cosine by 8
        winwid=8;
        winoff=4;
        weights = new double[winwid];
        dthet = 2*pi/winwid;
        for(int i=0; i<winwid; i++)  {
            weights[i]=(1-cos(i*dthet))/2;
        }
        setupTrans(winwid);
        break;
    case 2: // raised cos by 16
        winwid=16;
        winoff=6;
        weights = new double[winwid];
        dthet = 2*pi/(winwid-2);
        weights[0]=0.0;
        for(int i=1; i<winwid-1; i++)  {
            weights[i]=(1-cos((i-1)*dthet))/2;
        }
        weights[winwid-1]=0.0;
        //for(int i=0; i<winwid; i++)  weights[i]=1.0;
        setupTrans(winwid);
    case 3: // raised cos by 15
        winwid=15;
        winoff=15;
	setSize(15);
        weights = new double[winwid];
        dthet = 2*pi/(winwid-2);
        weights[0]=0.0;
        for(int i=1; i<winwid-1; i++)  {
            weights[i]=(1-cos((i-1)*dthet))/2;
        }
        weights[winwid-1]=0.0;
        for(int i=0; i<winwid; i++)  weights[i]=1.0;
        //for(int i=0; i<winwid; i++)  weights[i]=1.0;
        setupTrans(winwid);
    }
    weightBox->setCurrentIndex ( ibase );
}

void EnvWidget::setSize ( const size_t value )
{
    int add = value%winoff;
    size_t newval = value + add;
    if ( newval == N_ ) return;
    dirty = true;
    cStep = 0;
    if ( N_ != 0 ) {
        delete[] U_;
        delete[] X_;
        delete[] Ideal_;
        delete[] Usum;
        if(nbas) {
            delete[] Ub;
        }

        dirty = true;
    }
    N_ = newval;
    U_ = new double[N_];
    X_ = new double[N_];
    Usum = new double[N_];
    for(int i=0; i<N_; i++) Usum[i]=0.0;
    if(nbas) {
        Ub = new double[nbas*N_];
    }

    Ideal_ = new double[N_];
    initSin ( cycles );
}

void EnvWidget::step ( const size_t nStep )
{
    int istart;
    if( dirty ) {
        updateCoef(method);
        dirty=false;
    }
    if(totCFL == N_ / 2.0 ) {
        totCFL -= CFL;
        getIdeal();
        std::cout << "Initialize Ub = Ideal_ at " << totCFL << std::endl;
        int ustart=0;
        for( int iwin = winwid/2-1; iwin<N_ ; iwin+=winoff) {
            istart = iwin-(winwid)/2;
            if( istart < 0 ) istart= N_+istart;
            for( size_t i=0; i<winwid; i++) {
                gsl_vector_set(UVec,i,Ideal_[istart]*weights[i]);
                istart++;
                if(istart == N_) istart = 0;
            }
            // transform U  -  (U*window)-> Utran
            gsl_linalg_LU_solve(Mforw,permut,UVec,TranVec);
            for( size_t i=0; i<nbas; i++) {
                Ub[ustart] = gsl_vector_get(TranVec,i);
                ustart++;
            }
        }
        totCFL = N_/2.0;
    }
    if(unstable) return;
    samset = false;
    int ubstart=0;
    bool printinout=true;
    for ( size_t ns = 0; ns < nStep; ns++ ) {
        for( int iwin = 0; iwin<N_ ; iwin += winoff) {
            istart = iwin;
	    
            if( istart < 0 ) istart= N_+istart;
            int isw = istart+winwid;
            if(isw >= N_) isw -= N_;
            //gsl_vector_set(UVec,ipad,(U_[istart]+U_[isw])*weights[0]/2.0);
            for( size_t i=0; i<winwid; i++) {
                gsl_vector_set(UVec,i,U_[istart]*weights[i]);
                istart++;
                if(istart == N_) istart = 0;
            }
          
            // transform U  -  (U*window)-> Utran
            gsl_linalg_LU_solve(Mforw,permut,UVec,TranVec);
            // update Utran -
            // fillB
            gsl_blas_dgemv(CblasNoTrans,1.0,Right,TranVec,0.0,BVec);
            if(ntime == 3) {
                // make vector view from Ub
                UbView = gsl_vector_view_array (Ub+ubstart, nbas);
                // fill Cvec
                gsl_blas_dgemv(CblasNoTrans,1.0,FarRight,&UbView.vector,0.0,CVec);
                gsl_blas_daxpy(1.0,CVec,BVec);
                // save TranVec to Ub
                for( size_t i=0; i<nbas; i++) Ub[ubstart+i] = gsl_vector_get(TranVec,i);
                ubstart += nbas;
            }
            //    solve
            gsl_linalg_LU_solve(Left,lpermut,BVec,TranVec);
            // inverse transform TranVec
            gsl_blas_dgemv(CblasNoTrans,1.0,Mback,TranVec,0.0,UVec);
           
            // sum to new U
            istart = iwin;

            if( istart < 0 ) istart= N_+istart;
            for( size_t i=0; i<winwid; i++) {
                Usum[istart] += gsl_vector_get(UVec,i);
		if(Usum[istart] > 1e16) unstable = true;
                istart++;
                if(istart == N_) istart = 0;
            }

        }
        cStep++;
        totCFL += CFL;
        double * temp;
        temp = U_;
        U_ = Usum;
        Usum = temp;
        for(size_t n=0; n<N_; n++) Usum[n] = 0.0;
    }
}

void EnvWidget::setImplicit(double value) {
    if(value == impl) return;
    dirty = true;
    impl = value;
    beta = 1 - impl - back;
}

void EnvWidget::setBackward(double value) {
    if(value == back) return;
    dirty = true;
    back = value;
    beta = 1 -impl - back;
}

void EnvWidget::updateCoef(int value) {
    if( value != method ) setMethod(value);
    switch(method) {
    default:
        method = 0;
    case 0 : // LSP - central time
        lspc();
        break;
    case 1 : // LSP - Adams
        lspa();
        break;
    case 2 : // FEM - central time
        femc();
        break;
    case 3 : // FEM - Adams
        fema();
    }

}

void EnvWidget::setCFL(const double value) {
    if(CFL == value) return;
    dirty = true;
    CFL = value;
    //updateCoef(method);
}

void EnvWidget::setMethod(int index) {
    if( index == method ) return;
    dirty = true;
    method = index;
    ntime=2;

    switch( method ) {
    default:
        method = 0;
    case 0:
        setTitle( tr("Envelope - LSP"));
        if(back != 0.0) {
            impl = 0.5;
            beta = 0.5;
            back = 0.0;
        }
        break;
    case 1:
        setTitle( tr("Envelope - LSP - Adams"));
        ntime=3;
        if(back == 0.0) {
            impl = 5.0/12;
            beta = 2.0/3;
            back = -1.0/12;
        }
        break;
    case 2:
        setTitle( tr("Envelope - FEM"));
        if(back != 0.0) {
            impl = 0.5;
            beta = 0.5;
            back = 0.0;
        }
        break;
    case 3:
        setTitle( tr("Envelope - FEM - Adams"));
        ntime=3;
        if(back == 0.0) {
            impl = 5.0/12;
            beta = 2.0/3;
            back = -1.0/12;
        }
        break;
    }
    implInput->setValue(impl);
    backInput->setValue(back);
}

void EnvWidget::lspc() {
    double a0,a1,p2;
    //double den0,den1,den2,den3;
    //double cimp,cbet;

    a0 = impl * impl * CFL * CFL;//a2
    a1 = impl * beta * CFL * CFL;//ab
    //cimp = CFL*impl;
    //cbet = CFL*beta;
    p2 = pi * pi;
    //den0 = 1.0/(pi+2.0);
    //den1 = 1.0/(3*pi+16);
    //den2 = 1.0/(3*pi + 4);
    //den3 = 1.0/(3*pi + 8);
    switch(nbas) {
    default:
        allocate_gsl(4);
    case 4:
        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Left,0,1,0);
        gsl_matrix_set(Left,0,2,0);
        gsl_matrix_set(Left,0,3,0);

        gsl_matrix_set(Right,0,0,1);
        gsl_matrix_set(Right,0,1,0);
        gsl_matrix_set(Right,0,2,0);
        gsl_matrix_set(Right,0,3,0);


        gsl_matrix_set(Left,1,0,0);
        gsl_matrix_set(Left,1,1,0.25*p2*a0 + 1);
        gsl_matrix_set(Left,1,2,0);
        gsl_matrix_set(Left,1,3,0);

        gsl_matrix_set(Right,1,0,0);
        gsl_matrix_set(Right,1,1,-0.25*p2*a1 + 1);
        gsl_matrix_set(Right,1,2,0);
        gsl_matrix_set(Right,1,3,-0.5*pi*CFL);


        gsl_matrix_set(Left,2,0,0);
        gsl_matrix_set(Left,2,1,0);
        gsl_matrix_set(Left,2,2,p2*a0 + 1);
        gsl_matrix_set(Left,2,3,0);

        gsl_matrix_set(Right,2,0,0);
        gsl_matrix_set(Right,2,1,0);
        gsl_matrix_set(Right,2,2,-p2*a1 + 1);
        gsl_matrix_set(Right,2,3,0);


        gsl_matrix_set(Left,3,0,0);
        gsl_matrix_set(Left,3,1,0);
        gsl_matrix_set(Left,3,2,0);
        gsl_matrix_set(Left,3,3,0.25*p2*a0 + 1);

        gsl_matrix_set(Right,3,0,0);
        gsl_matrix_set(Right,3,1,0.5*pi*CFL);
        gsl_matrix_set(Right,3,2,0);
        gsl_matrix_set(Right,3,3,-0.25*p2*a1 + 1);
        break;
    case 8:
        gsl_matrix_set_all(Left,0.0);
        gsl_matrix_set_all(Right,0.0);
        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Right,0,0,1);
        gsl_matrix_set(Left,1,1,0.0625*p2*a0 + 1);
        gsl_matrix_set(Right,1,1,-0.0625*p2*a1 + 1);
        gsl_matrix_set(Right,1,5,-0.25*(pi)*CFL);
        gsl_matrix_set(Left,2,2,0.25*p2*a0 + 1);
        gsl_matrix_set(Right,2,2,-0.25*p2*a1 + 1);
        gsl_matrix_set(Right,2,6,-0.5*(pi)*CFL);
        gsl_matrix_set(Left,3,3,0.5625*p2*a0 + 1);
        gsl_matrix_set(Right,3,3,-0.5625*p2*a1 + 1);
        gsl_matrix_set(Right,3,7,-0.75*(pi)*CFL);
        gsl_matrix_set(Left,4,4,p2*a0 + 1);
        gsl_matrix_set(Right,4,4,-p2*a1 + 1);
        gsl_matrix_set(Left,5,5,0.0625*p2*a0 + 1);
        gsl_matrix_set(Right,5,1,0.25*(pi)*CFL);
        gsl_matrix_set(Right,5,5,-0.0625*p2*a1 + 1);
        gsl_matrix_set(Left,6,6,0.25*p2*a0 + 1);
        gsl_matrix_set(Right,6,2,0.5*(pi)*CFL);
        gsl_matrix_set(Right,6,6,-0.25*p2*a1 + 1);
        gsl_matrix_set(Left,7,7,0.5625*p2*a0 + 1);
        gsl_matrix_set(Right,7,3,0.75*(pi)*CFL);
        gsl_matrix_set(Right,7,7,-0.5625*p2*a1 + 1);
        break;
    case 15:
        gsl_matrix_set_all(Left,0.0);
        gsl_matrix_set_all(Right,0.0);


        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Right,0,0,1);
        gsl_matrix_set(Left,1,1,4.0/225*p2*a0 + 1);
        gsl_matrix_set(Right,1,1,-4.0/225*p2*a1 + 1);
        gsl_matrix_set(Right,1,8,-2.0/15*(pi)*CFL);
        gsl_matrix_set(Left,2,2,16.0/225*p2*a0 + 1);
        gsl_matrix_set(Right,2,2,-16.0/225*p2*a1 + 1);
        gsl_matrix_set(Right,2,9,-4.0/15*(pi)*CFL);
        gsl_matrix_set(Left,3,3,4.0/25*p2*a0 + 1);
        gsl_matrix_set(Right,3,3,-4.0/25*p2*a1 + 1);
        gsl_matrix_set(Right,3,10,-2.0/5*(pi)*CFL);
        gsl_matrix_set(Left,4,4,64.0/225*p2*a0 + 1);
        gsl_matrix_set(Right,4,4,-64.0/225*p2*a1 + 1);
        gsl_matrix_set(Right,4,11,-8.0/15*(pi)*CFL);
        gsl_matrix_set(Left,5,5,4.0/9*p2*a0 + 1);
        gsl_matrix_set(Right,5,5,-4.0/9*p2*a1 + 1);
        gsl_matrix_set(Right,5,12,-2.0/3*(pi)*CFL);
        gsl_matrix_set(Left,6,6,16.0/25*p2*a0 + 1);
        gsl_matrix_set(Right,6,6,-16.0/25*p2*a1 + 1);
        gsl_matrix_set(Right,6,13,-4.0/5*(pi)*CFL);
        gsl_matrix_set(Left,7,7,196.0/225*p2*a0 + 1);
        gsl_matrix_set(Right,7,7,-196.0/225*p2*a1 + 1);
        gsl_matrix_set(Right,7,14,-14.0/15*(pi)*CFL);
        gsl_matrix_set(Left,8,8,4.0/225*p2*a0 + 1);
        gsl_matrix_set(Right,8,1,2.0/15*(pi)*CFL);
        gsl_matrix_set(Right,8,8,-4.0/225*p2*a1 + 1);
        gsl_matrix_set(Left,9,9,16.0/225*p2*a0 + 1);
        gsl_matrix_set(Right,9,2,4.0/15*(pi)*CFL);
        gsl_matrix_set(Right,9,9,-16.0/225*p2*a1 + 1);
        gsl_matrix_set(Left,10,10,4.0/25*p2*a0 + 1);
        gsl_matrix_set(Right,10,3,2.0/5*(pi)*CFL);
        gsl_matrix_set(Right,10,10,-4.0/25*p2*a1 + 1);
        gsl_matrix_set(Left,11,11,64.0/225*p2*a0 + 1);
        gsl_matrix_set(Right,11,4,8.0/15*(pi)*CFL);
        gsl_matrix_set(Right,11,11,-64.0/225*p2*a1 + 1);
        gsl_matrix_set(Left,12,12,4.0/9*p2*a0 + 1);
        gsl_matrix_set(Right,12,5,2.0/3*(pi)*CFL);
        gsl_matrix_set(Right,12,12,-4.0/9*p2*a1 + 1);
        gsl_matrix_set(Left,13,13,16.0/25*p2*a0 + 1);
        gsl_matrix_set(Right,13,6,4.0/5*(pi)*CFL);
        gsl_matrix_set(Right,13,13,-16.0/25*p2*a1 + 1);
        gsl_matrix_set(Left,14,14,196.0/225*p2*a0 + 1);
        gsl_matrix_set(Right,14,7,14.0/15*(pi)*CFL);
        gsl_matrix_set(Right,14,14,-196.0/225*p2*a1 + 1);
        break;
    case 16:
        gsl_matrix_set_all(Left,0.0);
        gsl_matrix_set_all(Right,0.0);
        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Right,0,0,1);
        gsl_matrix_set(Left,1,1,0.015625*p2*a0 + 1);
        gsl_matrix_set(Right,1,1,-0.015625*p2*a1 + 1);
        gsl_matrix_set(Right,1,9,-0.125*(pi)*CFL);
        gsl_matrix_set(Left,2,2,0.0625*p2*a0 + 1);
        gsl_matrix_set(Right,2,2,-0.0625*p2*a1 + 1);
        gsl_matrix_set(Right,2,10,-0.25*(pi)*CFL);
        gsl_matrix_set(Left,3,3,0.140625*p2*a0 + 1);
        gsl_matrix_set(Right,3,3,-0.140625*p2*a1 + 1);
        gsl_matrix_set(Right,3,11,-0.375*(pi)*CFL);
        gsl_matrix_set(Left,4,4,0.25*p2*a0 + 1);
        gsl_matrix_set(Right,4,4,-0.25*p2*a1 + 1);
        gsl_matrix_set(Right,4,12,-0.5*(pi)*CFL);
        gsl_matrix_set(Left,5,5,25/64*p2*a0 + 1);
        gsl_matrix_set(Right,5,5,-25/64*p2*a1 + 1);
        gsl_matrix_set(Right,5,13,-0.625*(pi)*CFL);
        gsl_matrix_set(Left,6,6,9/16*p2*a0 + 1);
        gsl_matrix_set(Right,6,6,-9/16*p2*a1 + 1);
        gsl_matrix_set(Right,6,14,-0.25*(pi)*CFL);
        gsl_matrix_set(Left,7,7,49/64*p2*a0 + 1);
        gsl_matrix_set(Right,7,7,-49/64*p2*a1 + 1);
        gsl_matrix_set(Right,7,15,-0.875*(pi)*CFL);
        gsl_matrix_set(Left,8,8,p2*a0 + 1);
        gsl_matrix_set(Right,8,8,-p2*a1 + 1);
        gsl_matrix_set(Left,9,9,0.015625*p2*a0 + 1);
        gsl_matrix_set(Right,9,1,0.125*(pi)*CFL);
        gsl_matrix_set(Right,9,9,-0.015625*p2*a1 + 1);
        gsl_matrix_set(Left,10,10,0.0625*p2*a0 + 1);
        gsl_matrix_set(Right,10,2,0.25*(pi)*CFL);
        gsl_matrix_set(Right,10,10,-0.0625*p2*a1 + 1);
        gsl_matrix_set(Left,11,11,0.140625*p2*a0 + 1);
        gsl_matrix_set(Right,11,3,0.375*(pi)*CFL);
        gsl_matrix_set(Right,11,11,-0.140625*p2*a1 + 1);
        gsl_matrix_set(Left,12,12,0.25*p2*a0 + 1);
        gsl_matrix_set(Right,12,4,0.5*(pi)*CFL);
        gsl_matrix_set(Right,12,12,-0.25*p2*a1 + 1);
        gsl_matrix_set(Left,13,13,25/64*p2*a0 + 1);
        gsl_matrix_set(Right,13,5,0.625*(pi)*CFL);
        gsl_matrix_set(Right,13,13,-25/64*p2*a1 + 1);
        gsl_matrix_set(Left,14,14,9/16*p2*a0 + 1);
        gsl_matrix_set(Right,14,6,0.25*(pi)*CFL);
        gsl_matrix_set(Right,14,14,-9/16*p2*a1 + 1);
        gsl_matrix_set(Left,15,15,49/64*p2*a0 + 1);
        gsl_matrix_set(Right,15,7,0.875*(pi)*CFL);
        gsl_matrix_set(Right,15,15,-49/64*p2*a1 + 1);

    }


    gsl_matrix_set_all(FarRight,0.0);

    gsl_linalg_LU_decomp(Left,lpermut,&signum);

}


void EnvWidget::lspa() {
    double a0,a1,a2,p2,cib,cb;
    a0 = impl * impl * CFL * CFL;//a2
    a1 = impl * beta * CFL * CFL;//ab
    a2 = impl * back * CFL * CFL;//aB
    cib = CFL *( impl+beta );
    cb = CFL*back;
    p2 = pi * pi;

}


void EnvWidget::femc() {
    //double p2;
    //double den0,den1,den2,den3;
    double cimp,cbet;

    cimp = CFL*impl;
    cbet = CFL*beta;
    //p2 = pi * pi;
    //den0 = 1.0/(pi+2.0);
    //den1 = 1.0/(3*pi+16);
    //den2 = 1.0/(3*pi + 4);
    //den3 = 1.0/(3*pi + 8);
    switch(nbas) {
    default:
        allocate_gsl(4);
    case 4:
        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Left,0,1,0);
        gsl_matrix_set(Left,0,2,0);
        gsl_matrix_set(Left,0,3,0);

        gsl_matrix_set(Right,0,0,1);
        gsl_matrix_set(Right,0,1,0);
        gsl_matrix_set(Right,0,2,0);
        gsl_matrix_set(Right,0,3,0);


        gsl_matrix_set(Left,1,0,0);
        gsl_matrix_set(Left,1,1,1);
        gsl_matrix_set(Left,1,2,0);
        gsl_matrix_set(Left,1,3,0.5*pi*cimp);

        gsl_matrix_set(Right,1,0,0);
        gsl_matrix_set(Right,1,1,1);
        gsl_matrix_set(Right,1,2,0);
        gsl_matrix_set(Right,1,3,-0.5*pi*cbet);


        gsl_matrix_set(Left,2,0,0);
        gsl_matrix_set(Left,2,1,0);
        gsl_matrix_set(Left,2,2,1);
        gsl_matrix_set(Left,2,3,0);

        gsl_matrix_set(Right,2,0,0);
        gsl_matrix_set(Right,2,1,0);
        gsl_matrix_set(Right,2,2,1);
        gsl_matrix_set(Right,2,3,0);


        gsl_matrix_set(Left,3,0,0);
        gsl_matrix_set(Left,3,1,-0.5*pi*cimp);
        gsl_matrix_set(Left,3,2,0);
        gsl_matrix_set(Left,3,3,1);

        gsl_matrix_set(Right,3,0,0);
        gsl_matrix_set(Right,3,1,0.5*pi*cbet);
        gsl_matrix_set(Right,3,2,0);
        gsl_matrix_set(Right,3,3,1);
        break;
    case 8:
        gsl_matrix_set_all(Left,0.0);
        gsl_matrix_set_all(Right,0.0);
        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Right,0,0,1);
        gsl_matrix_set(Left,1,1,1);
        gsl_matrix_set(Left,1,5,0.25*pi*cimp);
        gsl_matrix_set(Right,1,1,1);
        gsl_matrix_set(Right,1,5,-0.25*pi*cbet);
        gsl_matrix_set(Left,2,2,1);
        gsl_matrix_set(Left,2,6,0.5*pi*cimp);
        gsl_matrix_set(Right,2,2,1);
        gsl_matrix_set(Right,2,6,-0.5*pi*cbet);
        gsl_matrix_set(Left,3,3,1);
        gsl_matrix_set(Left,3,7,0.75*pi*cimp);
        gsl_matrix_set(Right,3,3,1);
        gsl_matrix_set(Right,3,7,-0.75*pi*cbet);
        gsl_matrix_set(Left,4,4,1);
        gsl_matrix_set(Right,4,4,1);
        gsl_matrix_set(Left,5,1,-0.25*pi*cimp);
        gsl_matrix_set(Left,5,5,1);
        gsl_matrix_set(Right,5,1,0.25*pi*cbet);
        gsl_matrix_set(Right,5,5,1);
        gsl_matrix_set(Left,6,2,-0.5*pi*cimp);
        gsl_matrix_set(Left,6,6,1);
        gsl_matrix_set(Right,6,2,0.5*pi*cbet);
        gsl_matrix_set(Right,6,6,1);
        gsl_matrix_set(Left,7,3,-0.75*pi*cimp);
        gsl_matrix_set(Left,7,7,1);
        gsl_matrix_set(Right,7,3,0.75*pi*cbet);
        gsl_matrix_set(Right,7,7,1);
        break;
    case 15:
        gsl_matrix_set_all(Left,0.0);
        gsl_matrix_set_all(Right,0.0);


        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Right,0,0,1);
        gsl_matrix_set(Left,1,1,1);
        gsl_matrix_set(Left,1,8,2.0/15*pi*cimp);
        gsl_matrix_set(Right,1,1,1);
        gsl_matrix_set(Right,1,8,-2.0/15*pi*cbet);
        gsl_matrix_set(Left,2,2,1);
        gsl_matrix_set(Left,2,9,4.0/15*pi*cimp);
        gsl_matrix_set(Right,2,2,1);
        gsl_matrix_set(Right,2,9,-4.0/15*pi*cbet);
        gsl_matrix_set(Left,3,3,1);
        gsl_matrix_set(Left,3,10,2.0/5*pi*cimp);
        gsl_matrix_set(Right,3,3,1);
        gsl_matrix_set(Right,3,10,-2.0/5*pi*cbet);
        gsl_matrix_set(Left,4,4,1);
        gsl_matrix_set(Left,4,11,8.0/15*pi*cimp);
        gsl_matrix_set(Right,4,4,1);
        gsl_matrix_set(Right,4,11,-8.0/15*pi*cbet);
        gsl_matrix_set(Left,5,5,1);
        gsl_matrix_set(Left,5,12,2.0/3*pi*cimp);
        gsl_matrix_set(Right,5,5,1);
        gsl_matrix_set(Right,5,12,-2.0/3*pi*cbet);
        gsl_matrix_set(Left,6,6,1);
        gsl_matrix_set(Left,6,13,4.0/5*pi*cimp);
        gsl_matrix_set(Right,6,6,1);
        gsl_matrix_set(Right,6,13,-4.0/5*pi*cbet);
        gsl_matrix_set(Left,7,7,1);
        gsl_matrix_set(Left,7,14,14.0/15*pi*cimp);
        gsl_matrix_set(Right,7,7,1);
        gsl_matrix_set(Right,7,14,-14.0/15*pi*cbet);
        gsl_matrix_set(Left,8,1,-2.0/15*pi*cimp);
        gsl_matrix_set(Left,8,8,1);
        gsl_matrix_set(Right,8,1,2.0/15*pi*cbet);
        gsl_matrix_set(Right,8,8,1);
        gsl_matrix_set(Left,9,2,-4.0/15*pi*cimp);
        gsl_matrix_set(Left,9,9,1);
        gsl_matrix_set(Right,9,2,4.0/15*pi*cbet);
        gsl_matrix_set(Right,9,9,1);
        gsl_matrix_set(Left,10,3,-2.0/5*pi*cimp);
        gsl_matrix_set(Left,10,10,1);
        gsl_matrix_set(Right,10,3,2.0/5*pi*cbet);
        gsl_matrix_set(Right,10,10,1);
        gsl_matrix_set(Left,11,4,-8.0/15*pi*cimp);
        gsl_matrix_set(Left,11,11,1);
        gsl_matrix_set(Right,11,4,8.0/15*pi*cbet);
        gsl_matrix_set(Right,11,11,1);
        gsl_matrix_set(Left,12,5,-2.0/3*pi*cimp);
        gsl_matrix_set(Left,12,12,1);
        gsl_matrix_set(Right,12,5,2.0/3*pi*cbet);
        gsl_matrix_set(Right,12,12,1);
        gsl_matrix_set(Left,13,6,-4.0/5*pi*cimp);
        gsl_matrix_set(Left,13,13,1);
        gsl_matrix_set(Right,13,6,4.0/5*pi*cbet);
        gsl_matrix_set(Right,13,13,1);
        gsl_matrix_set(Left,14,7,-14.0/15*pi*cimp);
        gsl_matrix_set(Left,14,14,1);
        gsl_matrix_set(Right,14,7,14.0/15*pi*cbet);
        gsl_matrix_set(Right,14,14,1);
        break;
    case 16:
        gsl_matrix_set_all(Left,0.0);
        gsl_matrix_set_all(Right,0.0);
        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Right,0,0,1);
        gsl_matrix_set(Left,1,1,1);
        gsl_matrix_set(Left,1,9,0.125*pi*cimp);
        gsl_matrix_set(Right,1,1,1);
        gsl_matrix_set(Right,1,9,-0.125*pi*cbet);
        gsl_matrix_set(Left,2,2,1);
        gsl_matrix_set(Left,2,10,0.25*pi*cimp);
        gsl_matrix_set(Right,2,2,1);
        gsl_matrix_set(Right,2,10,-0.25*pi*cbet);
        gsl_matrix_set(Left,3,3,1);
        gsl_matrix_set(Left,3,11,0.375*pi*cimp);
        gsl_matrix_set(Right,3,3,1);
        gsl_matrix_set(Right,3,11,-0.375*pi*cbet);
        gsl_matrix_set(Left,4,4,1);
        gsl_matrix_set(Left,4,12,0.5*pi*cimp);
        gsl_matrix_set(Right,4,4,1);
        gsl_matrix_set(Right,4,12,-0.5*pi*cbet);
        gsl_matrix_set(Left,5,5,1);
        gsl_matrix_set(Left,5,13,0.625*pi*cimp);
        gsl_matrix_set(Right,5,5,1);
        gsl_matrix_set(Right,5,13,-0.625*pi*cbet);
        gsl_matrix_set(Left,6,6,1);
        gsl_matrix_set(Left,6,14,0.25*pi*cimp);
        gsl_matrix_set(Right,6,6,1);
        gsl_matrix_set(Right,6,14,-0.25*pi*cbet);
        gsl_matrix_set(Left,7,7,1);
        gsl_matrix_set(Left,7,15,0.875*pi*cimp);
        gsl_matrix_set(Right,7,7,1);
        gsl_matrix_set(Right,7,15,-0.875*pi*cbet);
        gsl_matrix_set(Left,8,8,1);
        gsl_matrix_set(Right,8,8,1);
        gsl_matrix_set(Left,9,1,-0.125*pi*cimp);
        gsl_matrix_set(Left,9,9,1);
        gsl_matrix_set(Right,9,1,0.125*pi*cbet);
        gsl_matrix_set(Right,9,9,1);
        gsl_matrix_set(Left,10,2,-0.25*pi*cimp);
        gsl_matrix_set(Left,10,10,1);
        gsl_matrix_set(Right,10,2,0.25*pi*cbet);
        gsl_matrix_set(Right,10,10,1);
        gsl_matrix_set(Left,11,3,-0.375*pi*cimp);
        gsl_matrix_set(Left,11,11,1);
        gsl_matrix_set(Right,11,3,0.375*pi*cbet);
        gsl_matrix_set(Right,11,11,1);
        gsl_matrix_set(Left,12,4,-0.5*pi*cimp);
        gsl_matrix_set(Left,12,12,1);
        gsl_matrix_set(Right,12,4,0.5*pi*cbet);
        gsl_matrix_set(Right,12,12,1);
        gsl_matrix_set(Left,13,5,-0.625*pi*cimp);
        gsl_matrix_set(Left,13,13,1);
        gsl_matrix_set(Right,13,5,0.625*pi*cbet);
        gsl_matrix_set(Right,13,13,1);
        gsl_matrix_set(Left,14,6,-0.25*pi*cimp);
        gsl_matrix_set(Left,14,14,1);
        gsl_matrix_set(Right,14,6,0.25*pi*cbet);
        gsl_matrix_set(Right,14,14,1);
        gsl_matrix_set(Left,15,7,-0.875*pi*cimp);
        gsl_matrix_set(Left,15,15,1);
        gsl_matrix_set(Right,15,7,0.875*pi*cbet);
        gsl_matrix_set(Right,15,15,1);

    }


    gsl_matrix_set_all(FarRight,0.0);

    gsl_linalg_LU_decomp(Left,lpermut,&signum);
}

void EnvWidget::fema() {
}

void EnvWidget::rk() {
    //double p2,den0,den1,den2,den3;

    //p2 = pi * pi;
    //den0 = 1.0/(pi+2.0);
    //den1 = 1.0/(3*pi+16);
    //den2 = 1.0/(3*pi + 4);
    //den3 = 1.0/(3*pi + 8);
    switch(nbas) {
    default:
        allocate_gsl(4);
    case 4:
        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Left,0,1,0);
        gsl_matrix_set(Left,0,2,0);
        gsl_matrix_set(Left,0,3,0);

        gsl_matrix_set(Right,0,0,0);
        gsl_matrix_set(Right,0,1,0);
        gsl_matrix_set(Right,0,2,0);
        gsl_matrix_set(Right,0,3,0);


        gsl_matrix_set(Left,1,0,0);
        gsl_matrix_set(Left,1,1,1);
        gsl_matrix_set(Left,1,2,0);
        gsl_matrix_set(Left,1,3,0);

        gsl_matrix_set(Right,1,0,0);
        gsl_matrix_set(Right,1,1,0);
        gsl_matrix_set(Right,1,2,0);
        gsl_matrix_set(Right,1,3,-0.5*pi*CFL);


        gsl_matrix_set(Left,2,0,0);
        gsl_matrix_set(Left,2,1,0);
        gsl_matrix_set(Left,2,2,1);
        gsl_matrix_set(Left,2,3,0);

        gsl_matrix_set(Right,2,0,0);
        gsl_matrix_set(Right,2,1,0);
        gsl_matrix_set(Right,2,2,0);
        gsl_matrix_set(Right,2,3,0);


        gsl_matrix_set(Left,3,0,0);
        gsl_matrix_set(Left,3,1,0);
        gsl_matrix_set(Left,3,2,0);
        gsl_matrix_set(Left,3,3,1);

        gsl_matrix_set(Right,3,0,0);
        gsl_matrix_set(Right,3,1,0.5*pi*CFL);
        gsl_matrix_set(Right,3,2,0);
        gsl_matrix_set(Right,3,3,0);
        break;
    case 8:
        gsl_matrix_set_all(Left,0.0);
        gsl_matrix_set_all(Right,0.0);
        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Left,1,1,1);
        gsl_matrix_set(Right,1,5,-0.25*pi*CFL);
        gsl_matrix_set(Left,2,2,1);
        gsl_matrix_set(Right,2,6,-0.5*pi*CFL);
        gsl_matrix_set(Left,3,3,1);
        gsl_matrix_set(Right,3,7,-0.75*pi*CFL);
        gsl_matrix_set(Left,4,4,1);
        gsl_matrix_set(Left,5,5,1);
        gsl_matrix_set(Right,5,1,0.25*pi*CFL);
        gsl_matrix_set(Left,6,6,1);
        gsl_matrix_set(Right,6,2,0.5*pi*CFL);
        gsl_matrix_set(Left,7,7,1);
        gsl_matrix_set(Right,7,3,0.75*pi*CFL);
        break;

    case 15:
        gsl_matrix_set_all(Left,0.0);
        gsl_matrix_set_all(Right,0.0);

        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Left,1,1,1);
        gsl_matrix_set(Right,1,8,-2.0/15*pi*CFL);
        gsl_matrix_set(Left,2,2,1);
        gsl_matrix_set(Right,2,9,-4.0/15*pi*CFL);
        gsl_matrix_set(Left,3,3,1);
        gsl_matrix_set(Right,3,10,-2.0/5*pi*CFL);
        gsl_matrix_set(Left,4,4,1);
        gsl_matrix_set(Right,4,11,-8.0/15*pi*CFL);
        gsl_matrix_set(Left,5,5,1);
        gsl_matrix_set(Right,5,12,-2.0/3*pi*CFL);
        gsl_matrix_set(Left,6,6,1);
        gsl_matrix_set(Right,6,13,-4.0/5*pi*CFL);
        gsl_matrix_set(Left,7,7,1);
        gsl_matrix_set(Right,7,14,-14.0/15*pi*CFL);
        gsl_matrix_set(Left,8,8,1);
        gsl_matrix_set(Right,8,1,2.0/15*pi*CFL);
        gsl_matrix_set(Left,9,9,1);
        gsl_matrix_set(Right,9,2,4.0/15*pi*CFL);
        gsl_matrix_set(Left,10,10,1);
        gsl_matrix_set(Right,10,3,2.0/5*pi*CFL);
        gsl_matrix_set(Left,11,11,1);
        gsl_matrix_set(Right,11,4,8.0/15*pi*CFL);
        gsl_matrix_set(Left,12,12,1);
        gsl_matrix_set(Right,12,5,2.0/3*pi*CFL);
        gsl_matrix_set(Left,13,13,1);
        gsl_matrix_set(Right,13,6,4.0/5*pi*CFL);
        gsl_matrix_set(Left,14,14,1);
        gsl_matrix_set(Right,14,7,14.0/15*pi*CFL);
	break;
    case 16:
        gsl_matrix_set_all(Left,0.0);
        gsl_matrix_set_all(Right,0.0);
        gsl_matrix_set(Left,0,0,1);
        gsl_matrix_set(Left,1,1,1);
        gsl_matrix_set(Right,1,9,-0.125*pi*CFL);
        gsl_matrix_set(Left,2,2,1);
        gsl_matrix_set(Right,2,10,-0.25*pi*CFL);
        gsl_matrix_set(Left,3,3,1);
        gsl_matrix_set(Right,3,11,-0.375*pi*CFL);
        gsl_matrix_set(Left,4,4,1);
        gsl_matrix_set(Right,4,12,-0.5*pi*CFL);
        gsl_matrix_set(Left,5,5,1);
        gsl_matrix_set(Right,5,13,-0.625*pi*CFL);
        gsl_matrix_set(Left,6,6,1);
        gsl_matrix_set(Right,6,14,-0.25*pi*CFL);
        gsl_matrix_set(Left,7,7,1);
        gsl_matrix_set(Right,7,15,-0.875*pi*CFL);
        gsl_matrix_set(Left,8,8,1);
        gsl_matrix_set(Left,9,9,1);
        gsl_matrix_set(Right,9,1,0.125*pi*CFL);
        gsl_matrix_set(Left,10,10,1);
        gsl_matrix_set(Right,10,2,0.25*pi*CFL);
        gsl_matrix_set(Left,11,11,1);
        gsl_matrix_set(Right,11,3,0.375*pi*CFL);
        gsl_matrix_set(Left,12,12,1);
        gsl_matrix_set(Right,12,4,0.5*pi*CFL);
        gsl_matrix_set(Left,13,13,1);
        gsl_matrix_set(Right,13,5,0.625*pi*CFL);
        gsl_matrix_set(Left,14,14,1);
        gsl_matrix_set(Right,14,6,0.25*pi*CFL);
        gsl_matrix_set(Left,15,15,1);
        gsl_matrix_set(Right,15,7,0.875*pi*CFL);
    }


    gsl_matrix_set_all(FarRight,0.0);

    gsl_linalg_LU_decomp(Left,lpermut,&signum);
}

void EnvWidget::initSin(const double value) {
    cStep = 0;
    unstable = false;
    if ( N_ == 0 ) setSize ( 100 );
    totCFL = N_ / 2.0;
    cycles = value;
    dx = 2 * pi / N_;
    for ( size_t i = 0; i < N_; i++ ) {
        X_[i] = dx * i;
        //U_[i] = sin ( cycles * X_[i] + travel);
    }
    getIdeal();
    for ( size_t i = 0; i < N_; i++ ) {
        U_[i] = Ideal_[i];
    }
}

double* EnvWidget::getU() {
    return U_;
}

void EnvWidget::allocate_gsl(int size) {
    if(size == nbas) return;
    if(nbas) {
        // deallocate matricies
        gsl_vector_free(UVec);
        gsl_vector_free(BVec);
        gsl_vector_free(CVec);
        gsl_vector_free(TranVec);
        gsl_matrix_free(Mforw);
        gsl_permutation_free(permut);
        gsl_matrix_free(Mback);
        gsl_matrix_free(Left);
        gsl_permutation_free(lpermut);
        gsl_matrix_free(Right);
        gsl_matrix_free(FarRight);
        if(N_) delete[] Ub;
    }
    nbas = size;
    UVec = gsl_vector_alloc(nbas);
    BVec = gsl_vector_alloc(nbas);
    CVec = gsl_vector_alloc(nbas);
    TranVec = gsl_vector_alloc(nbas);
    Mforw = gsl_matrix_alloc(nbas,nbas);
    permut = gsl_permutation_alloc (nbas);
    Mback =  gsl_matrix_alloc(nbas,nbas);
    Left = gsl_matrix_alloc(nbas,nbas);
    lpermut = gsl_permutation_alloc(nbas);
    Right = gsl_matrix_alloc(nbas,nbas);
    FarRight = gsl_matrix_alloc(nbas,nbas);
    if(N_) Ub= new double[N_*nbas];
}

void EnvWidget::setupTrans(int size) {
    int signum;
    //double rsq;
    //rsq = 1.0/sqrt(2.0);
    switch(size) {
    default:
        size=4;
    case 4:
        allocate_gsl(4);
        gsl_matrix_set(Mback,0,0,1);
        gsl_matrix_set(Mback,0,1,-1);
        gsl_matrix_set(Mback,0,2,1);
        gsl_matrix_set(Mback,0,3,0);

        gsl_matrix_set(Mback,1,0,1);
        gsl_matrix_set(Mback,1,1,0);
        gsl_matrix_set(Mback,1,2,-1);
        gsl_matrix_set(Mback,1,3,-1);

        gsl_matrix_set(Mback,2,0,1);
        gsl_matrix_set(Mback,2,1,1);
        gsl_matrix_set(Mback,2,2,1);
        gsl_matrix_set(Mback,2,3,0);

        gsl_matrix_set(Mback,3,0,1);
        gsl_matrix_set(Mback,3,1,0);
        gsl_matrix_set(Mback,3,2,-1);
        gsl_matrix_set(Mback,3,3,1);
        break;
    case 8:
        allocate_gsl(8);

        gsl_matrix_set(Mback,3,0,1.0);
        gsl_matrix_set(Mback,3,1,0.7071067811865475);
        gsl_matrix_set(Mback,3,2,0.0);
        gsl_matrix_set(Mback,3,3,-0.7071067811865475);
        gsl_matrix_set(Mback,3,4,-1.0);
        gsl_matrix_set(Mback,3,5,-0.7071067811865475);
        gsl_matrix_set(Mback,3,6,-1.0);
        gsl_matrix_set(Mback,3,7,-0.7071067811865475);

        gsl_matrix_set(Mback,4,0,1.0);
        gsl_matrix_set(Mback,4,1,1.0);
        gsl_matrix_set(Mback,4,2,1.0);
        gsl_matrix_set(Mback,4,3,1.0);
        gsl_matrix_set(Mback,4,4,1.0);
        gsl_matrix_set(Mback,4,5,0.0);
        gsl_matrix_set(Mback,4,6,0.0);
        gsl_matrix_set(Mback,4,7,0.0);

        gsl_matrix_set(Mback,5,0,1.0);
        gsl_matrix_set(Mback,5,1,0.7071067811865475);
        gsl_matrix_set(Mback,5,2,0.0);
        gsl_matrix_set(Mback,5,3,-0.7071067811865475);
        gsl_matrix_set(Mback,5,4,-1.0);
        gsl_matrix_set(Mback,5,5,0.7071067811865475);
        gsl_matrix_set(Mback,5,6,1.0);
        gsl_matrix_set(Mback,5,7,0.7071067811865475);

        gsl_matrix_set(Mback,6,0,1.0);
        gsl_matrix_set(Mback,6,1,0.0);
        gsl_matrix_set(Mback,6,2,-1.0);
        gsl_matrix_set(Mback,6,3,0.0);
        gsl_matrix_set(Mback,6,4,1.0);
        gsl_matrix_set(Mback,6,5,1.0);
        gsl_matrix_set(Mback,6,6,0.0);
        gsl_matrix_set(Mback,6,7,-1.0);

        gsl_matrix_set(Mback,7,0,1.0);
        gsl_matrix_set(Mback,7,1,-0.7071067811865475);
        gsl_matrix_set(Mback,7,2,0.0);
        gsl_matrix_set(Mback,7,3,0.7071067811865475);
        gsl_matrix_set(Mback,7,4,-1.0);
        gsl_matrix_set(Mback,7,5,0.7071067811865475);
        gsl_matrix_set(Mback,7,6,-1.0);
        gsl_matrix_set(Mback,7,7,0.7071067811865475);
        break;
    case 15:
        allocate_gsl(15);
        // location -7.5 = + 7.5
        gsl_matrix_set(Mback,0,0,1.0);
        gsl_matrix_set(Mback,0,1,-1.0);
        gsl_matrix_set(Mback,0,2,1.0);
        gsl_matrix_set(Mback,0,3,-1.0);
        gsl_matrix_set(Mback,0,4,1.0);
        gsl_matrix_set(Mback,0,5,-1.0);
        gsl_matrix_set(Mback,0,6,1.0);
        gsl_matrix_set(Mback,0,7,-1.0);
        gsl_matrix_set(Mback,0,8,0.0);
        gsl_matrix_set(Mback,0,9,2.288475490443933e-17);
        gsl_matrix_set(Mback,0,10,-3.432713235665899e-17);
        gsl_matrix_set(Mback,0,11,4.576950980887865e-17);
        gsl_matrix_set(Mback,0,12,-5.721188726109832e-17);
        gsl_matrix_set(Mback,0,13,6.865426471331798e-17);
        gsl_matrix_set(Mback,0,14,3.092566029697801e-17);

// location -6.5
        gsl_matrix_set(Mback,1,0,1.0);
        gsl_matrix_set(Mback,1,1,-0.9135454576426009);
        gsl_matrix_set(Mback,1,2,0.6691306063588583);
        gsl_matrix_set(Mback,1,3,-0.3090169943749476);
        gsl_matrix_set(Mback,1,4,-0.1045284632676532);
        gsl_matrix_set(Mback,1,5,0.5000000000000008);
        gsl_matrix_set(Mback,1,6,-0.8090169943749471);
        gsl_matrix_set(Mback,1,7,0.9781476007338058);
        gsl_matrix_set(Mback,1,8,-0.4067366430758001);
        gsl_matrix_set(Mback,1,9,0.7431448254773941);
        gsl_matrix_set(Mback,1,10,-0.9510565162951535);
        gsl_matrix_set(Mback,1,11,0.9945218953682734);
        gsl_matrix_set(Mback,1,12,-0.8660254037844382);
        gsl_matrix_set(Mback,1,13,0.5877852522924735);
        gsl_matrix_set(Mback,1,14,-0.2079116908177585);

// location -5.5
        gsl_matrix_set(Mback,2,0,1.0);
        gsl_matrix_set(Mback,2,1,-0.6691306063588581);
        gsl_matrix_set(Mback,2,2,-0.1045284632676538);
        gsl_matrix_set(Mback,2,3,0.8090169943749471);
        gsl_matrix_set(Mback,2,4,-0.9781476007338055);
        gsl_matrix_set(Mback,2,5,0.4999999999999996);
        gsl_matrix_set(Mback,2,6,0.3090169943749465);
        gsl_matrix_set(Mback,2,7,-0.9135454576426004);
        gsl_matrix_set(Mback,2,8,-0.7431448254773943);
        gsl_matrix_set(Mback,2,9,0.9945218953682733);
        gsl_matrix_set(Mback,2,10,-0.5877852522924735);
        gsl_matrix_set(Mback,2,11,-0.2079116908177600);
        gsl_matrix_set(Mback,2,12,0.8660254037844389);
        gsl_matrix_set(Mback,2,13,-0.9510565162951539);
        gsl_matrix_set(Mback,2,14,0.4067366430758014);

// location -4.5
        gsl_matrix_set(Mback,3,0,1.0);
        gsl_matrix_set(Mback,3,1,-0.3090169943749474);
        gsl_matrix_set(Mback,3,2,-0.8090169943749475);
        gsl_matrix_set(Mback,3,3,0.8090169943749475);
        gsl_matrix_set(Mback,3,4,0.3090169943749477);
        gsl_matrix_set(Mback,3,5,-1.0);
        gsl_matrix_set(Mback,3,6,0.3090169943749476);
        gsl_matrix_set(Mback,3,7,0.8090169943749471);
        gsl_matrix_set(Mback,3,8,-0.9510565162951536);
        gsl_matrix_set(Mback,3,9,0.5877852522924730);
        gsl_matrix_set(Mback,3,10,0.5877852522924730);
        gsl_matrix_set(Mback,3,11,-0.9510565162951535);
        gsl_matrix_set(Mback,3,12,-3.432713235665899e-17);
        gsl_matrix_set(Mback,3,13,0.9510565162951535);
        gsl_matrix_set(Mback,3,14,-0.5877852522924735);

// location -3.5
        gsl_matrix_set(Mback,4,0,1.0);
        gsl_matrix_set(Mback,4,1,0.1045284632676535);
        gsl_matrix_set(Mback,4,2,-0.9781476007338056);
        gsl_matrix_set(Mback,4,3,-0.3090169943749471);
        gsl_matrix_set(Mback,4,4,0.9135454576426009);
        gsl_matrix_set(Mback,4,5,0.5000000000000008);
        gsl_matrix_set(Mback,4,6,-0.8090169943749479);
        gsl_matrix_set(Mback,4,7,-0.6691306063588584);
        gsl_matrix_set(Mback,4,8,-0.9945218953682733);
        gsl_matrix_set(Mback,4,9,-0.2079116908177593);
        gsl_matrix_set(Mback,4,10,0.9510565162951537);
        gsl_matrix_set(Mback,4,11,0.4067366430758001);
        gsl_matrix_set(Mback,4,12,-0.8660254037844382);
        gsl_matrix_set(Mback,4,13,-0.5877852522924725);
        gsl_matrix_set(Mback,4,14,0.7431448254773941);

// location -2.5
        gsl_matrix_set(Mback,5,0,1.0);
        gsl_matrix_set(Mback,5,1,0.50);
        gsl_matrix_set(Mback,5,2,-0.4999999999999999);
        gsl_matrix_set(Mback,5,3,-1.0);
        gsl_matrix_set(Mback,5,4,-0.5000000000000002);
        gsl_matrix_set(Mback,5,5,0.4999999999999996);
        gsl_matrix_set(Mback,5,6,1.0);
        gsl_matrix_set(Mback,5,7,0.4999999999999996);
        gsl_matrix_set(Mback,5,8,-0.8660254037844386);
        gsl_matrix_set(Mback,5,9,-0.8660254037844387);
        gsl_matrix_set(Mback,5,10,0.0);
        gsl_matrix_set(Mback,5,11,0.8660254037844385);
        gsl_matrix_set(Mback,5,12,0.8660254037844389);
        gsl_matrix_set(Mback,5,13,2.288475490443933e-17);
        gsl_matrix_set(Mback,5,14,-0.8660254037844389);

// location -1.5
        gsl_matrix_set(Mback,6,0,1.0);
        gsl_matrix_set(Mback,6,1,0.8090169943749474);
        gsl_matrix_set(Mback,6,2,0.3090169943749474);
        gsl_matrix_set(Mback,6,3,-0.3090169943749477);
        gsl_matrix_set(Mback,6,4,-0.8090169943749475);
        gsl_matrix_set(Mback,6,5,-1.0);
        gsl_matrix_set(Mback,6,6,-0.8090169943749471);
        gsl_matrix_set(Mback,6,7,-0.3090169943749477);
        gsl_matrix_set(Mback,6,8,-0.5877852522924732);
        gsl_matrix_set(Mback,6,9,-0.9510565162951536);
        gsl_matrix_set(Mback,6,10,-0.9510565162951535);
        gsl_matrix_set(Mback,6,11,-0.5877852522924730);
        gsl_matrix_set(Mback,6,12,0.0);
        gsl_matrix_set(Mback,6,13,0.5877852522924736);
        gsl_matrix_set(Mback,6,14,0.9510565162951535);

// location -0.5
        gsl_matrix_set(Mback,7,0,1.0);
        gsl_matrix_set(Mback,7,1,0.9781476007338056);
        gsl_matrix_set(Mback,7,2,0.9135454576426009);
        gsl_matrix_set(Mback,7,3,0.8090169943749474);
        gsl_matrix_set(Mback,7,4,0.6691306063588582);
        gsl_matrix_set(Mback,7,5,0.5);
        gsl_matrix_set(Mback,7,6,0.3090169943749474);
        gsl_matrix_set(Mback,7,7,0.1045284632676535);
        gsl_matrix_set(Mback,7,8,-0.2079116908177593);
        gsl_matrix_set(Mback,7,9,-0.4067366430758002);
        gsl_matrix_set(Mback,7,10,-0.5877852522924732);
        gsl_matrix_set(Mback,7,11,-0.7431448254773942);
        gsl_matrix_set(Mback,7,12,-0.8660254037844386);
        gsl_matrix_set(Mback,7,13,-0.9510565162951536);
        gsl_matrix_set(Mback,7,14,-0.9945218953682733);

// location 0.5
        gsl_matrix_set(Mback,8,0,1.0);
        gsl_matrix_set(Mback,8,1,0.9781476007338056);
        gsl_matrix_set(Mback,8,2,0.9135454576426009);
        gsl_matrix_set(Mback,8,3,0.8090169943749474);
        gsl_matrix_set(Mback,8,4,0.6691306063588582);
        gsl_matrix_set(Mback,8,5,0.5);
        gsl_matrix_set(Mback,8,6,0.3090169943749474);
        gsl_matrix_set(Mback,8,7,0.1045284632676535);
        gsl_matrix_set(Mback,8,8,0.2079116908177593);
        gsl_matrix_set(Mback,8,9,0.4067366430758002);
        gsl_matrix_set(Mback,8,10,0.5877852522924732);
        gsl_matrix_set(Mback,8,11,0.7431448254773942);
        gsl_matrix_set(Mback,8,12,0.8660254037844386);
        gsl_matrix_set(Mback,8,13,0.9510565162951536);
        gsl_matrix_set(Mback,8,14,0.9945218953682733);

// location 1.5
        gsl_matrix_set(Mback,9,0,1.0);
        gsl_matrix_set(Mback,9,1,0.8090169943749474);
        gsl_matrix_set(Mback,9,2,0.3090169943749474);
        gsl_matrix_set(Mback,9,3,-0.3090169943749477);
        gsl_matrix_set(Mback,9,4,-0.8090169943749475);
        gsl_matrix_set(Mback,9,5,-1.0);
        gsl_matrix_set(Mback,9,6,-0.8090169943749471);
        gsl_matrix_set(Mback,9,7,-0.3090169943749477);
        gsl_matrix_set(Mback,9,8,0.5877852522924732);
        gsl_matrix_set(Mback,9,9,0.9510565162951536);
        gsl_matrix_set(Mback,9,10,0.9510565162951535);
        gsl_matrix_set(Mback,9,11,0.5877852522924730);
        gsl_matrix_set(Mback,9,12,0.0);
        gsl_matrix_set(Mback,9,13,-0.5877852522924736);
        gsl_matrix_set(Mback,9,14,-0.9510565162951535);

// location 2.5
        gsl_matrix_set(Mback,10,0,1.0);
        gsl_matrix_set(Mback,10,1,0.5);
        gsl_matrix_set(Mback,10,2,-0.4999999999999999);
        gsl_matrix_set(Mback,10,3,-1.0);
        gsl_matrix_set(Mback,10,4,-0.5000000000000002);
        gsl_matrix_set(Mback,10,5,0.4999999999999996);
        gsl_matrix_set(Mback,10,6,1.0);
        gsl_matrix_set(Mback,10,7,0.4999999999999996);
        gsl_matrix_set(Mback,10,8,0.8660254037844386);
        gsl_matrix_set(Mback,10,9,0.8660254037844387);
        gsl_matrix_set(Mback,10,10,0.0);
        gsl_matrix_set(Mback,10,11,-0.8660254037844385);
        gsl_matrix_set(Mback,10,12,-0.8660254037844389);
        gsl_matrix_set(Mback,10,13,-2.288475490443933e-17);
        gsl_matrix_set(Mback,10,14,0.8660254037844389);

// location 3.5
        gsl_matrix_set(Mback,11,0,1.0);
        gsl_matrix_set(Mback,11,1,0.1045284632676535);
        gsl_matrix_set(Mback,11,2,-0.9781476007338056);
        gsl_matrix_set(Mback,11,3,-0.3090169943749471);
        gsl_matrix_set(Mback,11,4,0.9135454576426009);
        gsl_matrix_set(Mback,11,5,0.5000000000000008);
        gsl_matrix_set(Mback,11,6,-0.8090169943749479);
        gsl_matrix_set(Mback,11,7,-0.6691306063588584);
        gsl_matrix_set(Mback,11,8,0.9945218953682733);
        gsl_matrix_set(Mback,11,9,0.2079116908177593);
        gsl_matrix_set(Mback,11,10,-0.9510565162951537);
        gsl_matrix_set(Mback,11,11,-0.4067366430758001);
        gsl_matrix_set(Mback,11,12,0.8660254037844382);
        gsl_matrix_set(Mback,11,13,0.5877852522924725);
        gsl_matrix_set(Mback,11,14,-0.7431448254773941);

// location 4.5
        gsl_matrix_set(Mback,12,0,1.0);
        gsl_matrix_set(Mback,12,1,-0.3090169943749474);
        gsl_matrix_set(Mback,12,2,-0.8090169943749475);
        gsl_matrix_set(Mback,12,3,0.8090169943749475);
        gsl_matrix_set(Mback,12,4,0.3090169943749477);
        gsl_matrix_set(Mback,12,5,-1.0);
        gsl_matrix_set(Mback,12,6,0.3090169943749476);
        gsl_matrix_set(Mback,12,7,0.8090169943749471);
        gsl_matrix_set(Mback,12,8,0.9510565162951536);
        gsl_matrix_set(Mback,12,9,-0.5877852522924730);
        gsl_matrix_set(Mback,12,10,-0.5877852522924730);
        gsl_matrix_set(Mback,12,11,0.9510565162951535);
        gsl_matrix_set(Mback,12,12,3.432713235665899e-17);
        gsl_matrix_set(Mback,12,13,-0.9510565162951535);
        gsl_matrix_set(Mback,12,14,0.5877852522924735);

// location 5.5
        gsl_matrix_set(Mback,13,0,1.0);
        gsl_matrix_set(Mback,13,1,-0.6691306063588581);
        gsl_matrix_set(Mback,13,2,-0.1045284632676538);
        gsl_matrix_set(Mback,13,3,0.8090169943749471);
        gsl_matrix_set(Mback,13,4,-0.9781476007338055);
        gsl_matrix_set(Mback,13,5,0.4999999999999996);
        gsl_matrix_set(Mback,13,6,0.3090169943749465);
        gsl_matrix_set(Mback,13,7,-0.9135454576426004);
        gsl_matrix_set(Mback,13,8,0.7431448254773943);
        gsl_matrix_set(Mback,13,9,-0.9945218953682733);
        gsl_matrix_set(Mback,13,10,0.5877852522924735);
        gsl_matrix_set(Mback,13,11,0.2079116908177600);
        gsl_matrix_set(Mback,13,12,-0.8660254037844389);
        gsl_matrix_set(Mback,13,13,0.9510565162951539);
        gsl_matrix_set(Mback,13,14,-0.4067366430758014);

// location 6.5
        gsl_matrix_set(Mback,14,0,1.0);
        gsl_matrix_set(Mback,14,1,-0.9135454576426009);
        gsl_matrix_set(Mback,14,2,0.6691306063588583);
        gsl_matrix_set(Mback,14,3,-0.3090169943749476);
        gsl_matrix_set(Mback,14,4,-0.1045284632676532);
        gsl_matrix_set(Mback,14,5,0.5000000000000008);
        gsl_matrix_set(Mback,14,6,-0.8090169943749471);
        gsl_matrix_set(Mback,14,7,0.9781476007338058);
        gsl_matrix_set(Mback,14,8,0.4067366430758001);
        gsl_matrix_set(Mback,14,9,-0.7431448254773941);
        gsl_matrix_set(Mback,14,10,0.9510565162951535);
        gsl_matrix_set(Mback,14,11,-0.9945218953682734);
        gsl_matrix_set(Mback,14,12,0.8660254037844382);
        gsl_matrix_set(Mback,14,13,-0.5877852522924735);
        gsl_matrix_set(Mback,14,14,0.2079116908177585);

    case 16:
        allocate_gsl(16);
        gsl_matrix_set(Mback,0,0,1.0);
        gsl_matrix_set(Mback,0,1,-1.0);
        gsl_matrix_set(Mback,0,2,1.0);
        gsl_matrix_set(Mback,0,3,-1.0);
        gsl_matrix_set(Mback,0,4,1.0);
        gsl_matrix_set(Mback,0,5,-1.0);
        gsl_matrix_set(Mback,0,6,1.0);
        gsl_matrix_set(Mback,0,7,-1.0);
        gsl_matrix_set(Mback,0,8,1.0);
        gsl_matrix_set(Mback,0,9,0.0);
        gsl_matrix_set(Mback,0,10,0.0);
        gsl_matrix_set(Mback,0,11,0.0);
        gsl_matrix_set(Mback,0,12,0.0);
        gsl_matrix_set(Mback,0,13,0.0);
        gsl_matrix_set(Mback,0,14,0.0);
        gsl_matrix_set(Mback,0,15,0.0);

        gsl_matrix_set(Mback,1,0,1.0);
        gsl_matrix_set(Mback,1,1,-0.9238795325112868);
        gsl_matrix_set(Mback,1,2,0.7071067811865475);
        gsl_matrix_set(Mback,1,3,-0.3826834323650898);
        gsl_matrix_set(Mback,1,4,0.0);
        gsl_matrix_set(Mback,1,5,0.3826834323650899);
        gsl_matrix_set(Mback,1,6,-0.7071067811865475);
        gsl_matrix_set(Mback,1,7,0.9238795325112868);
        gsl_matrix_set(Mback,1,8,-1.0);
        gsl_matrix_set(Mback,1,9,-0.3826834323650898);
        gsl_matrix_set(Mback,1,10,0.7071067811865475);
        gsl_matrix_set(Mback,1,11,-0.9238795325112868);
        gsl_matrix_set(Mback,1,12,1.0);
        gsl_matrix_set(Mback,1,13,-0.9238795325112867);
        gsl_matrix_set(Mback,1,14,0.7071067811865475);
        gsl_matrix_set(Mback,1,15,-0.3826834323650897);

        gsl_matrix_set(Mback,2,0,1.0);
        gsl_matrix_set(Mback,2,1,-0.7071067811865475);
        gsl_matrix_set(Mback,2,2,0.0);
        gsl_matrix_set(Mback,2,3,0.7071067811865475);
        gsl_matrix_set(Mback,2,4,-1.0);
        gsl_matrix_set(Mback,2,5,0.7071067811865475);
        gsl_matrix_set(Mback,2,6,0.0);
        gsl_matrix_set(Mback,2,7,-0.7071067811865475);
        gsl_matrix_set(Mback,2,8,1.0);
        gsl_matrix_set(Mback,2,9,-0.7071067811865475);
        gsl_matrix_set(Mback,2,10,1.0);
        gsl_matrix_set(Mback,2,11,-0.7071067811865475);
        gsl_matrix_set(Mback,2,12,0.0);
        gsl_matrix_set(Mback,2,13,0.7071067811865475);
        gsl_matrix_set(Mback,2,14,-1.0);
        gsl_matrix_set(Mback,2,15,0.7071067811865475);

        gsl_matrix_set(Mback,3,0,1.0);
        gsl_matrix_set(Mback,3,1,-0.3826834323650898);
        gsl_matrix_set(Mback,3,2,-0.7071067811865475);
        gsl_matrix_set(Mback,3,3,0.9238795325112868);
        gsl_matrix_set(Mback,3,4,0.0);
        gsl_matrix_set(Mback,3,5,-0.9238795325112868);
        gsl_matrix_set(Mback,3,6,0.7071067811865475);
        gsl_matrix_set(Mback,3,7,0.3826834323650899);
        gsl_matrix_set(Mback,3,8,-1.0);
        gsl_matrix_set(Mback,3,9,-0.9238795325112868);
        gsl_matrix_set(Mback,3,10,0.7071067811865475);
        gsl_matrix_set(Mback,3,11,0.3826834323650898);
        gsl_matrix_set(Mback,3,12,-1.0);
        gsl_matrix_set(Mback,3,13,0.3826834323650897);
        gsl_matrix_set(Mback,3,14,0.7071067811865475);
        gsl_matrix_set(Mback,3,15,-0.9238795325112867);

        gsl_matrix_set(Mback,4,0,1.0);
        gsl_matrix_set(Mback,4,1,0.0);
        gsl_matrix_set(Mback,4,2,-1.0);
        gsl_matrix_set(Mback,4,3,0.0);
        gsl_matrix_set(Mback,4,4,1.0);
        gsl_matrix_set(Mback,4,5,0.0);
        gsl_matrix_set(Mback,4,6,-1.0);
        gsl_matrix_set(Mback,4,7,0.0);
        gsl_matrix_set(Mback,4,8,1.0);
        gsl_matrix_set(Mback,4,9,-1.0);
        gsl_matrix_set(Mback,4,10,0.0);
        gsl_matrix_set(Mback,4,11,1.0);
        gsl_matrix_set(Mback,4,12,0.0);
        gsl_matrix_set(Mback,4,13,-1.0);
        gsl_matrix_set(Mback,4,14,0.0);
        gsl_matrix_set(Mback,4,15,1.0);

        gsl_matrix_set(Mback,5,0,1.0);
        gsl_matrix_set(Mback,5,1,0.3826834323650898);
        gsl_matrix_set(Mback,5,2,-0.7071067811865475);
        gsl_matrix_set(Mback,5,3,-0.9238795325112868);
        gsl_matrix_set(Mback,5,4,0.0);
        gsl_matrix_set(Mback,5,5,0.9238795325112868);
        gsl_matrix_set(Mback,5,6,0.7071067811865475);
        gsl_matrix_set(Mback,5,7,-0.3826834323650898);
        gsl_matrix_set(Mback,5,8,-1.0);
        gsl_matrix_set(Mback,5,9,-0.9238795325112868);
        gsl_matrix_set(Mback,5,10,-0.7071067811865475);
        gsl_matrix_set(Mback,5,11,0.3826834323650897);
        gsl_matrix_set(Mback,5,12,1.0);
        gsl_matrix_set(Mback,5,13,0.3826834323650898);
        gsl_matrix_set(Mback,5,14,-0.7071067811865475);
        gsl_matrix_set(Mback,5,15,-0.9238795325112868);

        gsl_matrix_set(Mback,6,0,1.0);
        gsl_matrix_set(Mback,6,1,0.7071067811865475);
        gsl_matrix_set(Mback,6,2,0.0);
        gsl_matrix_set(Mback,6,3,-0.7071067811865475);
        gsl_matrix_set(Mback,6,4,-1.0);
        gsl_matrix_set(Mback,6,5,-0.7071067811865475);
        gsl_matrix_set(Mback,6,6,0.0);
        gsl_matrix_set(Mback,6,7,0.7071067811865475);
        gsl_matrix_set(Mback,6,8,1.0);
        gsl_matrix_set(Mback,6,9,-0.7071067811865475);
        gsl_matrix_set(Mback,6,10,-1.0);
        gsl_matrix_set(Mback,6,11,-0.7071067811865475);
        gsl_matrix_set(Mback,6,12,0.0);
        gsl_matrix_set(Mback,6,13,0.7071067811865475);
        gsl_matrix_set(Mback,6,14,1.0);
        gsl_matrix_set(Mback,6,15,0.7071067811865475);

        gsl_matrix_set(Mback,7,0,1.0);
        gsl_matrix_set(Mback,7,1,0.9238795325112868);
        gsl_matrix_set(Mback,7,2,0.7071067811865475);
        gsl_matrix_set(Mback,7,3,0.3826834323650898);
        gsl_matrix_set(Mback,7,4,0.0);
        gsl_matrix_set(Mback,7,5,-0.3826834323650898);
        gsl_matrix_set(Mback,7,6,-0.7071067811865475);
        gsl_matrix_set(Mback,7,7,-0.9238795325112868);
        gsl_matrix_set(Mback,7,8,-1.0);
        gsl_matrix_set(Mback,7,9,-0.3826834323650898);
        gsl_matrix_set(Mback,7,10,-0.7071067811865475);
        gsl_matrix_set(Mback,7,11,-0.9238795325112868);
        gsl_matrix_set(Mback,7,12,-1.0);
        gsl_matrix_set(Mback,7,13,-0.9238795325112868);
        gsl_matrix_set(Mback,7,14,-0.7071067811865475);
        gsl_matrix_set(Mback,7,15,-0.3826834323650898);

        gsl_matrix_set(Mback,8,0,1.0);
        gsl_matrix_set(Mback,8,1,1.0);
        gsl_matrix_set(Mback,8,2,1.0);
        gsl_matrix_set(Mback,8,3,1.0);
        gsl_matrix_set(Mback,8,4,1.0);
        gsl_matrix_set(Mback,8,5,1.0);
        gsl_matrix_set(Mback,8,6,1.0);
        gsl_matrix_set(Mback,8,7,1.0);
        gsl_matrix_set(Mback,8,8,1.0);
        gsl_matrix_set(Mback,8,9,0.0);
        gsl_matrix_set(Mback,8,10,0.0);
        gsl_matrix_set(Mback,8,11,0.0);
        gsl_matrix_set(Mback,8,12,0.0);
        gsl_matrix_set(Mback,8,13,0.0);
        gsl_matrix_set(Mback,8,14,0.0);
        gsl_matrix_set(Mback,8,15,0.0);

        gsl_matrix_set(Mback,9,0,1.0);
        gsl_matrix_set(Mback,9,1,0.9238795325112868);
        gsl_matrix_set(Mback,9,2,0.7071067811865475);
        gsl_matrix_set(Mback,9,3,0.3826834323650898);
        gsl_matrix_set(Mback,9,4,0.0);
        gsl_matrix_set(Mback,9,5,-0.3826834323650898);
        gsl_matrix_set(Mback,9,6,-0.7071067811865475);
        gsl_matrix_set(Mback,9,7,-0.9238795325112868);
        gsl_matrix_set(Mback,9,8,-1.0);
        gsl_matrix_set(Mback,9,9,0.3826834323650898);
        gsl_matrix_set(Mback,9,10,0.7071067811865475);
        gsl_matrix_set(Mback,9,11,0.9238795325112868);
        gsl_matrix_set(Mback,9,12,1.0);
        gsl_matrix_set(Mback,9,13,0.9238795325112868);
        gsl_matrix_set(Mback,9,14,0.7071067811865475);
        gsl_matrix_set(Mback,9,15,0.3826834323650898);

        gsl_matrix_set(Mback,10,0,1.0);
        gsl_matrix_set(Mback,10,1,0.7071067811865475);
        gsl_matrix_set(Mback,10,2,0.0);
        gsl_matrix_set(Mback,10,3,-0.7071067811865475);
        gsl_matrix_set(Mback,10,4,-1.0);
        gsl_matrix_set(Mback,10,5,-0.7071067811865475);
        gsl_matrix_set(Mback,10,6,0.0);
        gsl_matrix_set(Mback,10,7,0.7071067811865475);
        gsl_matrix_set(Mback,10,8,1.0);
        gsl_matrix_set(Mback,10,9,0.7071067811865475);
        gsl_matrix_set(Mback,10,10,1.0);
        gsl_matrix_set(Mback,10,11,0.7071067811865475);
        gsl_matrix_set(Mback,10,12,0.0);
        gsl_matrix_set(Mback,10,13,-0.7071067811865475);
        gsl_matrix_set(Mback,10,14,-1.0);
        gsl_matrix_set(Mback,10,15,-0.7071067811865475);

        gsl_matrix_set(Mback,11,0,1.0);
        gsl_matrix_set(Mback,11,1,0.3826834323650898);
        gsl_matrix_set(Mback,11,2,-0.7071067811865475);
        gsl_matrix_set(Mback,11,3,-0.9238795325112868);
        gsl_matrix_set(Mback,11,4,0.0);
        gsl_matrix_set(Mback,11,5,0.9238795325112868);
        gsl_matrix_set(Mback,11,6,0.7071067811865475);
        gsl_matrix_set(Mback,11,7,-0.3826834323650898);
        gsl_matrix_set(Mback,11,8,-1.0);
        gsl_matrix_set(Mback,11,9,0.9238795325112868);
        gsl_matrix_set(Mback,11,10,0.7071067811865475);
        gsl_matrix_set(Mback,11,11,-0.3826834323650897);
        gsl_matrix_set(Mback,11,12,-1.0);
        gsl_matrix_set(Mback,11,13,-0.3826834323650898);
        gsl_matrix_set(Mback,11,14,0.7071067811865475);
        gsl_matrix_set(Mback,11,15,0.9238795325112868);

        gsl_matrix_set(Mback,12,0,1.0);
        gsl_matrix_set(Mback,12,1,0.0);
        gsl_matrix_set(Mback,12,2,-1.0);
        gsl_matrix_set(Mback,12,3,0.0);
        gsl_matrix_set(Mback,12,4,1.0);
        gsl_matrix_set(Mback,12,5,0.0);
        gsl_matrix_set(Mback,12,6,-1.0);
        gsl_matrix_set(Mback,12,7,0.0);
        gsl_matrix_set(Mback,12,8,1.0);
        gsl_matrix_set(Mback,12,9,1.0);
        gsl_matrix_set(Mback,12,10,0.0);
        gsl_matrix_set(Mback,12,11,-1.0);
        gsl_matrix_set(Mback,12,12,0.0);
        gsl_matrix_set(Mback,12,13,1.0);
        gsl_matrix_set(Mback,12,14,0.0);
        gsl_matrix_set(Mback,12,15,-1.0);

        gsl_matrix_set(Mback,13,0,1.0);
        gsl_matrix_set(Mback,13,1,-0.3826834323650898);
        gsl_matrix_set(Mback,13,2,-0.7071067811865475);
        gsl_matrix_set(Mback,13,3,0.9238795325112868);
        gsl_matrix_set(Mback,13,4,0.0);
        gsl_matrix_set(Mback,13,5,-0.9238795325112868);
        gsl_matrix_set(Mback,13,6,0.7071067811865475);
        gsl_matrix_set(Mback,13,7,0.3826834323650899);
        gsl_matrix_set(Mback,13,8,-1.0);
        gsl_matrix_set(Mback,13,9,0.9238795325112868);
        gsl_matrix_set(Mback,13,10,-0.7071067811865475);
        gsl_matrix_set(Mback,13,11,-0.3826834323650898);
        gsl_matrix_set(Mback,13,12,1.0);
        gsl_matrix_set(Mback,13,13,-0.3826834323650897);
        gsl_matrix_set(Mback,13,14,-0.7071067811865475);
        gsl_matrix_set(Mback,13,15,0.9238795325112867);

        gsl_matrix_set(Mback,14,0,1.0);
        gsl_matrix_set(Mback,14,1,-0.7071067811865475);
        gsl_matrix_set(Mback,14,2,0.0);
        gsl_matrix_set(Mback,14,3,0.7071067811865475);
        gsl_matrix_set(Mback,14,4,-1.0);
        gsl_matrix_set(Mback,14,5,0.7071067811865475);
        gsl_matrix_set(Mback,14,6,0.0);
        gsl_matrix_set(Mback,14,7,-0.7071067811865475);
        gsl_matrix_set(Mback,14,8,1.0);
        gsl_matrix_set(Mback,14,9,0.7071067811865475);
        gsl_matrix_set(Mback,14,10,-1.0);
        gsl_matrix_set(Mback,14,11,0.7071067811865475);
        gsl_matrix_set(Mback,14,12,0.0);
        gsl_matrix_set(Mback,14,13,-0.7071067811865475);
        gsl_matrix_set(Mback,14,14,1.0);
        gsl_matrix_set(Mback,14,15,-0.7071067811865475);

        gsl_matrix_set(Mback,15,0,1.0);
        gsl_matrix_set(Mback,15,1,-0.9238795325112868);
        gsl_matrix_set(Mback,15,2,0.7071067811865475);
        gsl_matrix_set(Mback,15,3,-0.3826834323650898);
        gsl_matrix_set(Mback,15,4,0.0);
        gsl_matrix_set(Mback,15,5,0.3826834323650899);
        gsl_matrix_set(Mback,15,6,-0.7071067811865475);
        gsl_matrix_set(Mback,15,7,0.9238795325112868);
        gsl_matrix_set(Mback,15,8,-1.0);
        gsl_matrix_set(Mback,15,9,0.3826834323650898);
        gsl_matrix_set(Mback,15,10,-0.7071067811865475);
        gsl_matrix_set(Mback,15,11,0.9238795325112868);
        gsl_matrix_set(Mback,15,12,-1.0);
        gsl_matrix_set(Mback,15,13,0.9238795325112867);
        gsl_matrix_set(Mback,15,14,-0.7071067811865475);
        gsl_matrix_set(Mback,15,15,0.3826834323650897);
        break;


    }
    gsl_matrix_memcpy(Mforw,Mback);
    gsl_linalg_LU_decomp(Mforw,permut,&signum);
}
bool EnvWidget::canSolve(int equ) {
    return pdeSolvers[equ][method];
}

