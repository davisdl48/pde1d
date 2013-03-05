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

#include "specwidget.h"
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <cmath>

SpecWidget::SpecWidget ( QWidget* parent ) : SolvWidget(parent)
{   int iwid = 4;
    setTitle(tr("Spectral FFT Method"));

    methodLabel = new QLabel ( tr ( "Method" ) );
    methodBox = new QComboBox ( this );
    methodBox->addItem( tr("LSP - Central time"));
    methodBox->addItem( tr("LSP - Adams"));
    methodBox->addItem( tr("FEM - Central time"));
    methodBox->addItem( tr("FEM - Adams"));
    methodBox->addItem( tr("RK-Mid Pt"));
    methodBox->addItem( tr("RK-Trap"));
    methodBox->addItem( tr("RK4"));
    methodBox->addItem( tr("Arb-RK4"));
    methodBox->addItem( tr("Phaser"));
    connect ( methodBox, SIGNAL( activated(int) ), this, SLOT( setMethod(int) ) );
    verticalLayout->insertWidget ( iwid++, methodLabel );
    verticalLayout->insertWidget ( iwid++, methodBox );
    /*
    weightLabel = new QLabel ( tr ( "Basis Functions" ) );
    weightBox = new QComboBox ( this );
    weightBox->addItem ( tr ( "Linear by 4" ) );
    weightBox->addItem ( tr ( "Cosine by 8" ) );
    weightBox->addItem ( tr ( "Cosine by 16" ) );
    weightBox->addItem ( tr ( "Cosine by 15" ) );
    verticalLayout->insertWidget ( iwid++, weightLabel );
    verticalLayout->insertWidget ( iwid++, weightBox );
    connect ( weightBox, SIGNAL( activated(int) ), this, SLOT( setBasis(int) ) );
    */
    plotNameEdit->setText ( title );
    implLabel = new QLabel ( tr ( "Implicit" ) );
    implInput = new KDoubleNumInput ( 0.0, 1.0, 5/12.0, parent, 1e-14, 14 );
    implInput->setValue(5/12.0);
    connect ( implInput, SIGNAL ( valueChanged ( double ) ), this, SLOT ( setImplicit ( double ) ) );
    verticalLayout->insertWidget ( iwid++, implLabel );
    verticalLayout->insertWidget ( iwid++, implInput );
    backLabel = new QLabel ( tr ( "Backward" ) );
    backInput = new KDoubleNumInput ( -1.0, 1.0, -1/12.0, parent, 1e-14, 14 );
    backInput->setValue(-1/12.0);
    connect ( backInput, SIGNAL ( valueChanged ( double ) ), this, SLOT ( setBackward(double)) );
    verticalLayout->insertWidget ( iwid++, backLabel );
    verticalLayout->insertWidget ( iwid++, backInput );

    /* Set the default values */
    setColor ( Qt::magenta );
    narrays=0;
    nStage = 0;
    method=-1;
    setMethod(2);
    setImplicit(0.5);
    //set_default_options(&options);
}

SpecWidget::SpecWidget ( const SpecWidget& other )
{

}

SpecWidget::~SpecWidget()
{
    if ( N_ != 0 ) {
        gsl_fft_real_wavetable_free (real_g);
        gsl_fft_halfcomplex_wavetable_free (hc_g);
        gsl_fft_real_workspace_free (work_g);
    }
    freeData();
}

SpecWidget& SpecWidget::operator= ( const SpecWidget& other )
{
    return *this;
}

bool SpecWidget::operator== ( const SpecWidget& other ) const
{
    return (this == &other);
}


void SpecWidget::setSize ( const size_t value )
{
    if ( value == N_ ) return;
    cStep = 0;
    if ( N_ != 0 ) {
        delete[] U_;
        delete[] X_;
        delete[] Ideal_;
        gsl_fft_real_wavetable_free (real_g);
        gsl_fft_halfcomplex_wavetable_free (hc_g);
        gsl_fft_real_workspace_free (work_g);

    }
    N_ = value;
    U_ = new double[N_];
    X_ = new double[N_];

    work_g = gsl_fft_real_workspace_alloc (N_);
    real_g = gsl_fft_real_wavetable_alloc (N_);
    hc_g = gsl_fft_halfcomplex_wavetable_alloc (N_);

    Ideal_ = new double[N_];
    initSin ( cycles );
    allocateData(nStage);
}

void SpecWidget::step ( const size_t nStep )
{
    double temp;
    if(totCFL == N_ / 2.0 && ntime == 3) {
        totCFL -= CFL;
        getIdeal();
        std::cout << "Initialize Ub = Ideal_ at " << totCFL << std::endl;
        for( size_t nn=0; nn<N_; nn++) data[2][nn] = Ideal_[nn];
        // transform data[2]
        totCFL = N_/2.0;
    }
    for ( size_t ns = 0; ns < nStep; ns++ ) {
        // copy U_ to data[1]
        for( size_t nn=0; nn<N_; nn++) data[1][nn] = U_[nn];


        // solve for data[0] -- ie new U_
        data[0][0] = data[1][0]; // always no updata to ~U[0]
        switch( method ) {
        default:
        case 0:
            lspc();
            break;
        case 1:
            lspa();
            break;
        case 2:
            femc();
            break;
        case 3:
            fema();
            break;
        case 4:
        case 5:
        case 6:
        case 7:
            for(int stg=0; stg < nStage; stg++) {
                for(size_t nn=0; nn<N_; nn++) data[1][nn] = U_[nn];
                for(int s=0 ; s<stg; s++) {
                    temp=b_a[nStage*stg+s];
                    if(temp != 0.0) for (size_t nn=0; nn<N_; nn++) data[1][nn] += temp*data[ ntime+s ][nn];
                }
                rk(data[ ntime+stg ],data[1]);
            }

            for(size_t nn = 0; nn <  N_ ; nn++ ) data[0][nn] = U_[nn];
            for(int j=0 ; j<nStage; j++) {
                temp = b_b[j];
                if( temp != 0 ) {
                    for(size_t nn = 0; nn <  N_ ; nn++ ) {
                        data[0][nn] += temp*data[ ntime+j ][nn];
                    }
                }
            }
            break;

        case 8:
            phaser();
        }


        // Update U_
        for( size_t nn=0; nn<N_; nn++) U_[nn] = data[0][nn];

        cStep++;
        totCFL += CFL;
    }
}

void SpecWidget::setImplicit(double value) {
    if(value == impl) return;
    impl = value;
    beta = 1 - impl - back;
}

void SpecWidget::setBackward(double value) {
    if(value == back) return;
    back = value;
    beta = 1 -impl - back;
}

void SpecWidget::setCFL(const double value) {
    if(CFL == value) return;
    CFL = value;
    //updateCoef(method);
}

void SpecWidget::setMethod(int index) {
    if( index == method ) return;
    method = index;
    ntime=2;

    switch( method ) {
    default:
        method = 0;
    case 0:
    case 1:
        setTitle( tr("Spectral - LSP"));
        if(back != 0.0) {
            impl = 0.5;
            beta = 0.5;
            back = 0.0;
        }
        break;
        setTitle( tr("Spectral - LSP - Adams"));
        ntime=3;
        if(back == 0.0) {
            impl = 5.0/12;
            beta = 2.0/3;
            back = -1.0/12;
        }
        break;
    case 2:
    case 3:
        setTitle( tr("Spectral - FEM"));
        if(back != 0.0) {
            impl = 0.5;
            beta = 0.5;
            back = 0.0;
        }
        break;
        setTitle( tr("Spectral - FEM - Adams"));
        ntime=3;
        if(back == 0.0) {
            impl = 5.0/12;
            beta = 2.0/3;
            back = -1.0/12;
        }
        break;
    case 4:// Mid-Point
//	y1 = y0 + h f( y0+h/2 f(y0))
//	nStage = 2
//	b_a = matrix(SR,[[0,0],[1/2,0]])
//	b_b=[0,1]
//	b_c=[0,1/2]
        setTitle( tr("Spectral -Mid-Point"));
        setNStage(2);
        b_a[nStage]=0.5;
        b_b[1]=1.0;
        b_c[1]=0.5;
        break;
    case 5: // Trapezoid
        //y1=y0+h/2f(y0) + h/2 f(y0+h f(y0))
//	nStage = 2
//	b_a = matrix(SR,[[0,0,],[1,0])
//	b_b = [1/2,1/2]
//	b_c=[0,1]
        setTitle( tr("Spectral -Trapezoid"));
        setNStage(2);
        b_a[nStage]=1.0;
        b_b[0]=0.5;
        b_b[1]=0.5;
        b_c[1]=1.0;

        break;
    case 6: // RK4
        /*
        * k_i = f(Y_i)
        * Y_i = y_0 + h sum a_ij k_j
        * y_1 = y_0 + h sum b_i k_i
        *
        b_a = matrix(SR,[[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]])
        b_b = vector(QQ,[1/6,1/3,1/3,1/6])
        b_c = vector(QQ,[0,1/2,1/2,1])
        b_k = vector(SR,[0,0,0,0])
        nStage = 4

        */
        setTitle( tr("Spectral -RK4"));
        setNStage(4);
        b_a[4]=0.5;
        b_a[9]=0.5;
        b_a[14]=1;
        b_b[0]=1.0/6.0;
        b_b[1]=1.0/3.0;
        b_b[2]=1.0/3.0;
        b_b[3]=1.0/6.0;
        b_c[1]=0.5;
        b_c[2]=0.5;
        b_c[3]=1;
        break;
    case 7: // Arb-RK-4
        /*
        * k_i = f(Y_i)
        * Y_i = y_0 + h sum a_ij k_j
        * y_1 = y_0 + h sum b_i k_i
        *
        b_a = matrix(SR,[[0,0,0,0],[1/2,0,0,0],[0,1/2,0,0],[0,0,1,0]])
        b_b = vector(QQ,[1/6,1/3,1/3,1/6])
        b_c = vector(QQ,[0,1/2,1/2,1])
        b_k = vector(SR,[0,0,0,0])
        nStage = 4

        */
        setTitle( tr("Spectral -Arb-RK-4"));
        setNStage(4);
        b_a[4]=0.25;
        b_a[9]=1.0/3.0;
        b_a[14]=0.5;
        b_b[3]=1.0;
        b_c[1]=0.25;
        b_c[2]=1.0/3.0;
        b_c[3]=0.5;
        break;
    case 8:
        setTitle( tr("Phaser"));
        break;
    }
    implInput->setValue(impl);
    backInput->setValue(back);
    methodBox->setCurrentIndex(method);
}

void SpecWidget::lspc() {
    double aa,ab;
    double dscal,ds2;
    double det;
    //double den0,den1,den2,den3;
    //double cimp,cbet;
    double ui;
    double ur;
    double rur,rui;
    dscal = 2*pi/N_;

    ds2=dscal*dscal;
    double right[2][2];
    double left[2][2];
    int np;

    aa = impl * impl * CFL * CFL;//a2
    ab = impl * beta * CFL * CFL;//ab
    // transform data[1]
    gsl_fft_real_transform (data[1], 1, N_, real_g, work_g);


    data[0][0] = data[1][0];
    for(size_t nn=0; nn<N_/2; nn++) {
        ur = data[1][2*nn+1];
        if( 2*nn+2 == N_) {
            data[0][2*nn+1] = ur*(1-ds2*aa*nn*nn)/(1+ds2*aa*nn*nn);
            //data[0][2*nn-1] = ur;
            continue;
        }
        ui = data[1][2*nn+2];
        np = nn+1;

        right[0][0] = 1-ds2*aa*np*np;
        right[0][1] = dscal*np;
        right[1][0] = -dscal*np;
        right[1][1] = 1-ds2*aa*np*np;
        rur=right[0][0]*ur+right[0][1]*ui;
        rui=right[1][0]*ur+right[1][1]*ui;

        left[0][0] = 1+ds2*aa*np*np;
        left[0][1] = 0.0;
        left[1][0] = 0.0;
        left[1][1] = 1+ds2*aa*np*np;
        det = left[0][0]*left[1][1]-left[0][1]*left[1][0];
        data[0][2*nn+1] = (left[1][1]*rur-left[0][1]*rui)/det;
        data[0][2*nn+2] = (-left[1][0]*rur+left[0][0]*rui)/det;
    }

    gsl_fft_halfcomplex_inverse (data[0], 1, N_, hc_g, work_g);
}


void SpecWidget::lspa() {
    double a0,a1,a2,p2,cib,cb;
    a0 = impl * impl * CFL * CFL;//a2
    a1 = impl * beta * CFL * CFL;//ab
    a2 = impl * back * CFL * CFL;//aB
    cib = CFL *( impl+beta );
    cb = CFL*back;
    p2 = pi * pi;      // transform data[1]
    gsl_fft_real_transform (data[1], 1, N_, real_g, work_g);
    for(size_t nn=0; nn<N_; nn++) data[0][nn] = data[1][nn];

    gsl_fft_halfcomplex_inverse (data[0], 1, N_, hc_g, work_g);

}


void SpecWidget::femc() {
    double dscal;
    double det;
    double cimp,cbet;
    //double den0,den1,den2,den3;
    //double cimp,cbet;
    double ui,ur;
    double rui,rur;

    double right[2][2];
    double left[2][2];
    dscal = 2*pi/N_;

    cimp = CFL*impl;
    cbet = CFL*beta;
    // transform data[1]
    gsl_fft_real_transform (data[1], 1, N_, real_g, work_g);


    data[0][0] = data[1][0];
    for(size_t nn=0; nn<N_/2; nn++) {
        ur = data[1][2*nn+1];
        if( 2*nn+2 == N_) {
            data[0][2*nn+1] = ur;
            continue;
        }
        ui = data[1][2*nn+2];

        right[0][0] = 1;
        right[0][1] = dscal*cbet*(nn+1);
        right[1][0] = -dscal*cbet*(nn+1);
        right[1][1] = 1;

        left[0][0] = 1;
        left[0][1] = -dscal*cimp*(nn+1);
        left[1][0] = dscal*cimp*(nn+1);
        left[1][1] = 1;

        det = left[0][0]*left[1][1]-left[0][1]*left[1][0];
        rur=right[0][0]*ur+right[0][1]*ui;
        rui=right[1][0]*ur+right[1][1]*ui;
        data[0][2*nn+1] = (left[1][1]*rur-left[0][1]*rui)/det;
        data[0][2*nn+2] = (-left[1][0]*rur+left[0][0]*rui)/det;
    }

    gsl_fft_halfcomplex_inverse (data[0], 1, N_, hc_g, work_g);
}

void SpecWidget::fema() {
    // transform data[1]
    gsl_fft_real_transform (data[1], 1, N_, real_g, work_g);
    for(size_t nn=0; nn<N_; nn++) data[0][nn] = data[1][nn];

    gsl_fft_halfcomplex_inverse (data[0], 1, N_, hc_g, work_g);

}

void SpecWidget::rk(double* kout,double * uin) {
    double dscal;
    double ui,ur;
    int np;

    dscal = 2*pi/N_;
    // transform data[1]
    gsl_fft_real_transform (uin, 1, N_, real_g, work_g);

    kout[0] = 0.0;
    for(size_t nn=0; nn<N_/2; nn++) {
        np=nn+1;
        ur = uin[2*nn+1];
        if( 2*nn+2 == N_) {
            kout[2*nn+1] = 0.0;
            continue;
        }
        ui = uin[2*nn+2];

        kout[2*nn+1] = ui*dscal*CFL*np;
        kout[2*nn+2] = -ur*dscal*CFL*np;
    }
    gsl_fft_halfcomplex_inverse (kout, 1, N_, hc_g, work_g);
}

void SpecWidget::initSin(const double value) {
    cStep = 0;
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

double* SpecWidget::getU() {
    return U_;
}

void SpecWidget::allocateData(size_t rksize) { // implicit ntime and N_ usage
    size_t newsize;
    if(N_ == 0) return;
    newsize = ntime+rksize;
    if(newsize != narrays) { //deallocate all data arrays
        freeData();
        data = new double*[newsize];
        narrays=newsize;
    }
    if(dsize != N_) { //reallocate all data arrays
        if(dsize) {
            for( size_t n =0; n<narrays; n++ ) {
                delete[] data[n];
            }
        }
        for( size_t n =0; n<narrays; n++ ) {
            data[n] = new double[N_];
            if( data[n] == 0 ) {
                std::cout << "data allocation failed - please close\n";
                N_ = 0;
                return;
            }
        }
        dsize = N_;
    }
}

void SpecWidget::freeData() {
    if(N_ == 0) return;
    for( size_t n =0; n<narrays; n++ ) {
        delete[] data[n];
    }
    dsize = 0;
    if(narrays) delete[] data;
    narrays=0;
}

void SpecWidget::setNStage(int arg1) {
    if(nStage != arg1) {
        if(nStage != 0) {
            delete[] b_a;
            delete[] b_b;
            delete[] b_c;
        }
        b_a = new double[arg1*arg1];
        b_b = new double[arg1];
        b_c = new double[arg1];
        if( N_ != 0) allocateData(arg1);
        nStage = arg1;
    }
    for(int i=0; i<nStage; i++ ) {
        for(int j=0; j<nStage; j++) b_a[j*nStage+i] = 0.0;
        b_b[i] = 0.0;
        b_c[i] = 0.0;
    }
}

void SpecWidget::phaser() {
    double dscal;
    double ui,ur;
    double dr,di;
    double dot;
    double mag,ph;
    int np;
    dscal = 2*pi/N_;
    // transform data[1]
    gsl_fft_real_transform (data[1], 1, N_, real_g, work_g);

    data[0][0] = data[1][0];
    for(size_t nn=0; nn<N_/2; nn++) {
        np = nn+1;
        ur = data[1][2*nn+1];
        if( 2*nn+2 == N_) {
            data[0][2*nn+1] = 0;
            continue;
        }
        ui = data[1][2*nn+2];
        mag = hypot(ui,ur);
        ph = atan2(ui,ur);

        dr = ui*dscal*CFL*np;
        di = -ur*dscal*CFL*np;

        if(mag > 1e-12) {
            dot= (dr*ur+di*ui)/mag;

            dr -= ur*dot/mag;
            di -= ui*dot/mag;
            if(ur*di-ui*dr > 0.0) {
                /// ??? positive rotation
                ph += hypot(dr,di)/mag;
            } else {
                ph -= hypot(dr,di)/mag;
            }
            mag += dot;

            ur = mag*cos(ph);
            ui = mag*sin(ph);
        } else {
            ur += dr;
            ui += di;
        }
        data[0][2*nn+1] = ur;
        data[0][2*nn+2] = ui;
    }
    gsl_fft_halfcomplex_inverse (data[0], 1, N_, hc_g, work_g);
}

