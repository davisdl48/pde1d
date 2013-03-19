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

#include "impwidget.h"
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <cmath>

ImpWidget::ImpWidget ( QWidget* parent ) : SolvWidget(parent)
{   int iwid = 4;
    setTitle(tr("LSP Central Time"));

    methodLabel = new QLabel ( tr ( "Method" ) );
    methodBox = new QComboBox ( this );
    methodBox->addItem( tr("LSP - Central time"));
    methodBox->addItem( tr("LSP - Adams"));
    methodBox->addItem( tr("FEM - Central time"));
    methodBox->addItem( tr("FEM - Adams"));
    connect ( methodBox, SIGNAL( activated(int) ), this, SLOT( setMethod(int) ) );
    verticalLayout->insertWidget ( iwid++, methodLabel );
    verticalLayout->insertWidget ( iwid++, methodBox );
    weightLabel = new QLabel ( tr ( "Basis Functions" ) );
    weightBox = new QComboBox ( this );
    weightBox->addItem ( tr ( "Linear Basis" ) );
    weightBox->addItem ( tr ( "Cosine Basis" ) );
    weightBox->addItem ( tr ( "Parabolic Basis" ) );
    weightBox->addItem ( tr ( "Linear Fourier" ) );
    verticalLayout->insertWidget ( iwid++, weightLabel );
    verticalLayout->insertWidget ( iwid++, weightBox );
    connect ( weightBox, SIGNAL( activated(int) ), this, SLOT( setBasis(int) ) );
    plotNameEdit->setText ( title );
    implLabel = new QLabel ( tr ( "Implicit" ) );
    implInput = new MyDoubInput ( 5/12.0, this, 0.0, 1.0,  1e-14, 14 );
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

    /* Set the default values for options argument: */
    options.Fact = DOFACT;
    options.Equil = YES;
    options.ColPerm = COLAMD;
    options.DiagPivotThresh = 1.0;
    options.Trans = NOTRANS;
    options.IterRefine = NOREFINE;
    options.SymmetricMode = NO;
    options.PivotGrowth = NO;
    options.ConditionNumber = NO;
    options.PrintStat = NO;
    /* Set the default values */
    setColor ( Qt::magenta );
    lwork = 0;
    nrhs  = 1;
    ncoef = 0;
    aexist = false;
    dirty = true;
    nblock = 1;
    matSize_gsl=0;
    transform = false;
    ibase=-1;
    method=-1;
    setBasis(0);
    setMethod(0);
  unstable = false;
    //set_default_options(&options);
}

ImpWidget::ImpWidget ( const ImpWidget& other )
{

}

ImpWidget::~ImpWidget()
{
    if ( N_ != 0 ) {

        delete[] Ub;
        delete[] Utran;

        SUPERLU_FREE (rhsb);
        SUPERLU_FREE (rhsx);
        SUPERLU_FREE (etree);
        SUPERLU_FREE (perm_r);
        SUPERLU_FREE (perm_c);
        SUPERLU_FREE (R);
        SUPERLU_FREE (C);
        SUPERLU_FREE (ferr);
        SUPERLU_FREE (berr);
        /// ???
        if(aexist) {
            /// ??? Destroy_CompCol_Matrix(&A);
            delete[] a;
            delete[] xa;
            delete[] asub;
        }
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        if ( lwork == 0 && !dirty) {
            Destroy_SuperNode_Matrix(&L);
            Destroy_CompCol_Matrix(&Up);
        } else if ( lwork > 0 ) {
            SUPERLU_FREE(work);
        }
    }

    if( matSize_gsl ) {
        // deallocate matricies
        gsl_vector_free(UVec);
        gsl_vector_free(TranVec);
        gsl_matrix_free(Mforw);
        gsl_matrix_free(Mback);
    }
}

ImpWidget& ImpWidget::operator= ( const ImpWidget& other )
{
    return *this;
}

bool ImpWidget::operator== ( const ImpWidget& other ) const
{
    return (this == &other);
}


void ImpWidget::setBasis ( int index )
{
    if(ibase != index) dirty = true;
    ibase = index;
    switch(ibase) {
    default:
        ibase = 0;
    case 0: // linear
    case 1: // raised cosine
        nvar = 3;
        nblock = 1;
        ivar = 1;
        break;
    case 2: // peicewise parabolic
        nvar = 5;
        nblock = 2;
        ivar=1;
        break;
    case 3: // Fourier 4
        nvar = 12;
        nblock = 4;
        ivar = 4;
        setupTrans();
        break;
    }
    ndn = ivar;
    nup = nvar - ivar;
    weightBox->setCurrentIndex ( ibase );
}

void ImpWidget::setSize ( const size_t value )
{
    int add = value%nblock;
    size_t newval = value + add;
    if ( newval == N_ ) return;
    dirty = true;
    cStep = 0;
    if ( N_ != 0 ) {
        delete[] U_;
        delete[] Ub;
        delete[] X_;
        delete[] Ideal_;
        delete[] Utran;

        SUPERLU_FREE (rhsb);
        SUPERLU_FREE (rhsx);
        SUPERLU_FREE (etree);
        SUPERLU_FREE (perm_r);
        SUPERLU_FREE (perm_c);
        SUPERLU_FREE (R);
        SUPERLU_FREE (C);
        SUPERLU_FREE (ferr);
        SUPERLU_FREE (berr);
        /// ???
        if(aexist) {
            /// ??? Destroy_CompCol_Matrix(&A);
            //delete[] a;
            //delete[] xa;
            //delete[] asub;
        }
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        if ( lwork == 0 && !dirty) {
            Destroy_SuperNode_Matrix(&L);
            Destroy_CompCol_Matrix(&Up);
        } else if ( lwork > 0 ) {
            SUPERLU_FREE(work);
        }
        aexist= false;
        dirty = true;
    }
    N_ = newval;
    U_ = new double[N_];
    Ub = new double[N_];
    X_ = new double[N_];
    Utran = new double[N_];

    Ideal_ = new double[N_];
    initSin ( cycles );

    if ( !(rhsb = doubleMalloc(N_)) ) ABORT("Malloc fails for rhsb[].");
    if ( !(rhsx = doubleMalloc(N_)) ) ABORT("Malloc fails for rhsx[].");
    dCreate_Dense_Matrix(&B, N_, 1, rhsb, N_, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&X, N_, 1, rhsx, N_, SLU_DN, SLU_D, SLU_GE);

    if ( !(etree = intMalloc(N_)) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_r = intMalloc(N_)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(N_)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(R = (double *) SUPERLU_MALLOC(N_ * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (double *) SUPERLU_MALLOC(N_ * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (double *) SUPERLU_MALLOC( sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (double *) SUPERLU_MALLOC( sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for berr[].");
}

void ImpWidget::step ( const size_t nStep )
{
    DNformat *Bstore;
    DNformat *Xstore;
    if(totCFL == N_ / 2.0 ) {
        totCFL += -CFL;
        getIdeal();
        std::cout << "Initialize Ub = Ideal_ at " << totCFL << std::endl;
        if(transform) {
            forwardTrans(Ideal_,Ub);
            forwardTrans(U_,Utran);
        } else {
            for(size_t n = 0; n < N_; n++) {
                Ub[n]=Ideal_[n];
            }
        }
        totCFL = N_/2.0;
    }
    if(unstable) return;
    if( dirty ) {
        if(aexist) {
            Destroy_CompCol_Matrix(&A);
            if ( lwork == 0 ) {
                Destroy_SuperNode_Matrix(&L);
                Destroy_CompCol_Matrix(&Up);
            } else if ( lwork > 0 ) {
                SUPERLU_FREE(work);
            }
            // these may be freed in dgssvx or Destroy_CompCol_Matrix I think
            aexist = false;
        }
        a = new double[nvar*N_];
        xa = new int[N_+1];
        asub = new int[nvar*N_];
        updateCoef(method);
        if(nblock == 1) {
            // load with coef[] ???
            fillA();
            std::cout << " nnz = " << nnz << std::endl;
        } else if (nup+ndn == 0) {
            // solve seperate independent blocks???
        } else {
            // load block system
            blockFillA();
        }

        dCreate_CompCol_Matrix(&A, N_, N_, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
        aexist = true;
        /* Initialize the statistics variables. */
        StatInit(&stat);

        dPrint_CompCol_Matrix("A matrix", &A);
        options.Fact = DOFACT;
        //options.ColPerm=NATURAL;
        options.ColPerm=COLAMD;
        options.PivotGrowth = NO;
        options.ConditionNumber = NO;
        /* ONLY PERFORM THE LU DECOMPOSITION */
        B.ncol = 0;  /* Indicate not to solve the system */
        dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
               &L, &Up, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
               &mem_usage, &stat, &info);
        //dPrint_CompCol_Matrix("A matrix", &A);
        printf("LU factorization: dgssvx() returns info %d\n", info);
        if ( info == 0 || info == N_+1 ) {

            if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", rpg);
            if ( options.ConditionNumber )
                printf("Recip. condition number = %e\n", rcond);

            printf("L\\U_ MB %.3f\ttotal MB needed %.3f\n",
                   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
            fflush(stdout);
            options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
            B.ncol = 1;
            dirty = false;
        } else if ( info > 0 && lwork == -1 ) {
            printf("** Estimated memory: %d bytes\n", info - n);
        }
        if ( options.PrintStat ) StatPrint(&stat);
        StatFree(&stat);
    }
    for ( size_t n = 0; n < nStep; n++ ) {
        ///set B matrix
        Bstore= (DNformat *)B.Store;
        rhsb=(double*)Bstore->nzval;
        if (nup+ndn == 0) {
            // solve seperate independent blocks???
            std::cout << "Discontinuous Galerkin Not Implimented\n";
            return;
        } else {
            // solve block system
            fillB();
        }

        //solve factored system
        StatInit(&stat);
        options.Fact = FACTORED;
        dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
               &L, &Up, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
               &mem_usage, &stat, &info);

        if( n == 1 ) printf("Triangular solve: dgssvx() returns info %d\n", info);

        //update U_
        if ( info == 0 || info == N_+1 ) {

            /* This is how you could access the solution matrix. */
            Xstore = (DNformat *)X.Store;
            rhsx = (double*)(Xstore->nzval);
            if(transform) {
                for ( size_t i = 0; i <  N_ ; i++ ) {
                    Ub[i] = Utran[i];
                    Utran[i]=rhsx[i];
		    if(Utran[i] > 1e16) unstable = true;
                }
            } else {
                for ( size_t i = 0; i <  N_ ; i++ ) {
                    Ub[i] = U_[i];
                    U_[i]=rhsx[i];
		    if(U_[i] > 1e16) unstable = true;
                }
            }
        } else {
            std::cout << "ERROR: Matrix Solution Failed   info = " << info << std::endl;
        }
        StatFree(&stat);
        cStep++;
        totCFL += CFL;
    }
}

void ImpWidget::setImplicit(double value) {
    if(value == impl) return;
    dirty = true;
    impl = value;
    beta = 1 - impl - back;
    implInput->setValue(impl);
}

void ImpWidget::setBackward(double value) {
    if(value == back) return;
    dirty = true;
    back = value;
    beta = 1 -impl - back;
    backInput->setValue(back);
}

void ImpWidget::updateCoef(int value) {
    if( value != method ) setMethod(value);
    if(ncoef != nblock*nvar*ntime) {
        if( ncoef != 0 ) {
            delete[] coef;
        }
        ncoef = nblock*nvar*ntime;
        coef = new double[ncoef];
    }
    switch(method) {
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
    std::cout << " Coefficients\n";
    for(int i= 0; i<ncoef; i++) {
        std::cout << "  "<< coef[i] ;
    }
    std::cout << std::endl;

}

void ImpWidget::setCFL(const double value) {
    if(CFL == value) return;
    dirty = true;
    CFL = value;
    //updateCoef(method);
}

void ImpWidget::fillA() {
    /*
    * nvar = 5
    * nblock = 1
    * ivar = 2
    * N_ = 8
    * nup = (nvar-ivar) = 3
    * ndn = ivar = 2
    *
    *   *0*     *1*     *2*     *3*     *4*     *5*     *6*     *7*
    *                  nup-1                           N_-ndn    N_-1
    *
    *  c[2]    c[3]    c[4]      0       0       0     c[0]    c[1]    *0*
    *  c[1]    c[2]    c[3]    c[4]      0       0       0     c[0]    *1*  ndn-1
    *  c[0]    c[1]    c[2]    c[3]    c[4]      0       0       0     *2*
    *    0     c[0]    c[1]    c[2]    c[3]    c[4]      0       0     *3*
    *    0       0     c[0]    c[1]    c[2]    c[3]    c[4]      0     *4*
    *    0       0       0     c[0]    c[1]    c[2]    c[3]    c[4]    *5*  N_-nup
    *  c[4]      0       0       0     c[0]    c[1]    c[2]    c[3]    *6*
    *  c[3]    c[4]      0       0       0     c[0]    c[1]    c[2]    *7*
    *
    *         nup-2
    *
    */
    int col,row;
    int cbeg;
    double temp;
    xa[0] = 0;
    nnz=0;
    //std::cout << " --col-  -row-  -cbeg-\n";
    for(col = 0; col < N_; col++) {
        if(col < nup ) { // upper left
            row = 0;
            cbeg = ivar + col;
        } else if( col >= N_-ndn ) { // upper right
            row = 0;
            cbeg = col -(N_ - ndn);
            //std::cout << col << "  " << row  << "  " << cbeg << std::endl;
        } else { // middle
            row = col - nup +1;
            cbeg = nvar-1;
        }
        //std::cout << col << "  " << row << "  "<< nvar << "  " << cbeg << std::endl;
        for(; cbeg >= 0 ; cbeg--) {
            temp = coef[cbeg];
            if(temp != 0.0 ) {
                a[nnz]=temp;
                asub[nnz] = row;
                nnz ++;
            }
            row++;
            if(row >= N_) break;
        }
        cbeg = 0;
        if(col < (nup-1)) { // lower left
            row = col + N_ - nup +1 ;
            cbeg = nvar - 1;
        } else if( col >= N_-ndn ) { // lower right
            row = col - nup +1;
            cbeg = nvar-1;
        }
        if(cbeg) {
            std::cout << col << "  " << row  << "  " << cbeg << std::endl;
            for(; cbeg >= 0 ; cbeg--) {
                temp = coef[cbeg];
                if(temp != 0.0 ) {
                    a[nnz]=temp;
                    asub[nnz] = row;
                    nnz ++;
                }
                row++;
                if(row >= N_) break;
            }
        }
        xa[col+1] = nnz;
    }
}

void ImpWidget::blockFillA() {
    /*
    * nvar = 5
    * nblock = 2
    * ivar = 2
    * N_ = 8
    * ndn = ivar
    * nup = nvar-ivar
    *
    *  *0*     *1*  |   *2*     *3*  |   *4*     *5*   |  *6*     *7*
    *                  nup-1                             N_-ndn    N_-1
    *
    * c0[2]   c0[3] |  c0[4]      0  |     0       0   |  c0[0]   c0[1]    *0*
    * c1[2]   c1[3] |  c1[4]      0  |     0       0   |  c1[0]   c1[1]    *1*
    * ---------------------------------------------------------------------
    * c0[0]   c0[1] |  c0[2]   c0[3] |  c0[4]      0   |     0       0     *2*
    * c1[0]   c1[1] |  c1[2]   c1[3] |  c1[4]      0   |     0       0     *3*
    * ---------------------------------------------------------------------
    *    0       0  |  c0[0]   c0[1] |  c0[2]   c0[3]  |  c0[4]      0     *4*
    *    0       0  |  c1[0]   c1[1] |  c1[2]   c1[3]  |  c1[4]      0     *5*
    * ---------------------------------------------------------------------
    * c0[4]      0  |     0       0  |  c0[0]   c0[1]  |  c0[2]   c0[3]    *6*
    * c1[4]      0  |     0       0  |  c1[0]   c1[1]  |  c1[2]   c1[3]    *7*
    *
    */
    int bnup = nup/nblock;
    int bndn = (ndn+1)/nblock;
    int bivar = (ivar+nblock-1)/nblock;
    int bN = N_/nblock;
    int bcol,bbeg;// block column, block row, block begin
    int bcolcol; // column within block;
    int bnvar=bndn + bnup;
    int col,row;
    int ccol;  // coef column within block
    double temp;

    xa[0] = 0;
    nnz=0;
    for(col = 0; col < N_; col++) {
        bcol = col/nblock;
        bcolcol = col - bcol*nblock;
        if(bcol < bnup ) { // upper left
            row = 0;
            bbeg = bivar + bcol;
        } else if( bcol >= bN-bndn ) { // upper right
            row = 0;
            bbeg = bcol -(bN - bndn);
        } else { // middle
            row = (bcol - bnup +1)*nblock;
            bbeg = bnvar-1;
        }
        //std::cout << col << "  " << row << "  "<< nvar << "  " << bbeg << std::endl;
        for(; bbeg >= 0 ; bbeg--) {
            ccol = (bbeg-bivar)*nblock + bcolcol + ivar;
            if(ccol < 0) continue;
            if(ccol >= nvar) continue;
            for(int nb=0; nb< nblock; nb++ ) {
                temp = coef[ccol+nb*nvar*ntime];
                if(temp != 0.0 ) {
                    a[nnz]=temp;
                    asub[nnz] = row;
                    nnz ++;
                }
                row++;
                if(row >= N_) break;
            }
            if(row >= N_) break;
        }
        bbeg = 0;
        if(bcol < (bnup-1)) { // lower left
            row = (bcol + bN - bnup +1) *nblock;
            bbeg = bnvar - 1;
        } else if( bcol >= bN-bndn ) { // lower right
            row = (bcol - bnup +1)*nblock;
            bbeg = bnvar-1;
        }
        if(bbeg) {
            for(; bbeg >= 0 ; bbeg--) {
                ccol = (bbeg-bivar)*nblock + bcolcol + ivar;
                if(ccol < 0) break;
                if(ccol >= nvar) continue;
                for(int nb=0; nb< nblock; nb++ ) {
                    temp = coef[ccol+nb*nvar*ntime];
                    if(temp != 0.0 ) {
                        a[nnz]=temp;
                        asub[nnz] = row;
                        nnz ++;
                    }
                    row++;
                    if(row >= N_) break;
                }
                if(row >= N_) break;
            }
        }
        xa[col+1] = nnz;
    }
}

void ImpWidget::fillB() {
    int ui;
    if(ibase == 3) { // Fourier Bases
        //need to translate U_ to Utran for each block
    }
    for( int i =0; i<N_; i += nblock ) {
        for(int nb = 0; nb < nblock; nb++) {
            rhsb[i+nb] = 0;
            ui = i-ndn;
            if(ui < 0) ui = N_+ui;
            //std::cout << i << " U_[n]  \n";
            for( int j =0; j<nvar; j++ ) {
                //std::cout << ui << "  ";
                if(transform) {
                    rhsb[i+nb] += coef[nb*nvar*ntime + nvar + j]*Utran[ui];
                } else {
                    rhsb[i+nb] += coef[nb*nvar*ntime + nvar + j]*U_[ui];
                }
                //std::cout << "coef [ "<< nb*nvar*ntime + nvar + j << " ] * U_[ " << ui << " ]  = ";
                //std::cout << coef[nb*nvar*ntime + nvar + j]<< " * " << U_[ui] << std::endl;
                ui++;
                if(ui == N_) ui=0;
            }
            //for(int nb = 0; nb < nblock; nb++) {
            //std::cout << "b[ " << i+nb << " ] = " << rhsb[i+nb] << std::endl;
            //}
        }
    }
    if(ntime == 2) return;
    for( int i =0; i<N_; i+=nblock ) {
        for(int nb = 0; nb < nblock; nb++) {
            ui = i-ndn;
            if(ui < 0) ui = N_+ui;
            for( int j =0; j<nvar; j++) {
                rhsb[i+nb] += coef[nb*nvar*ntime + 2*nvar + j]*Ub[ui];
            }
            ui++;
            if(ui == N_) ui=0;
        }
    }
}


void ImpWidget::setMethod(int index) {
    if( index == method ) return;
    transform = false;
    dirty = true;
    method = index;
    ntime=2;
    switch( method ) {
    default:
        method = 0;
    case 0:
        setTitle( tr("LSP - Central"));
        if(back != 0.0) {
            impl = 0.5;
            beta = 0.5;
            back = 0.0;
        }
        break;
    case 1:
        setTitle( tr("LSP - Adams"));
        ntime=3;
        if(back == 0.0) {
            impl = 5.0/12;
            beta = 2.0/3;
            back = -1.0/12;
        }
        break;
    case 2:
        setTitle( tr("FEM - Central"));
        if(back != 0.0) {
            impl = 0.5;
            beta = 0.5;
            back = 0.0;
        }
        break;
    case 3:
        setTitle( tr("FEM - Adams"));
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

void ImpWidget::lspc() {

    double a0,a1,p2;

    a0 = impl * impl * CFL * CFL;//a2
    a1 = impl * beta * CFL * CFL;//ab
    p2 = pi * pi;
    transform = false;
    switch(ibase) {
    default:
    case 0: // LSP - central time - linear
        /*[-C^2*a^2 + 1/6, 2*C^2*a^2 + 2/3, -C^2*a^2 + 1/6, -C^2*a*b - 1/2*C*a -
        1/2*C*b - 1/6, 2*C^2*a*b - 2/3, -C^2*a*b + 1/2*C*a + 1/2*C*b - 1/6]*/
        coef[0]= -a0 + 1.0/6.0;
        coef[1]=  2.0*a0 + 2.0/3.0;
        coef[2]=  -a0 + 1.0/6.0;

        coef[3]=  a1 + 0.5*CFL + 1.0/6.0;
        coef[4]=  -2.0*a1 + 2.0/3.0;
        coef[5]=  a1 - 0.5*CFL + 1.0/6.0;
        break;
    case 1: // LSP - central time - raised cosine
        /*[-1/8*pi^2*C^2*a^2 + 1/8, 1/4*pi^2*C^2*a^2 + 3/4, -1/8*pi^2*C^2*a^2 +
        1/8, -1/8*pi^2*C^2*a*b - 1/2*C*a - 1/2*C*b - 1/8, 1/4*pi^2*C^2*a*b -
        3/4, -1/8*pi^2*C^2*a*b + 1/2*C*a + 1/2*C*b - 1/8] */

        coef[0]= -0.125*p2*a0 + 0.125;
        coef[1]=  0.25*p2*a0 + 0.75;
        coef[2]=  -0.125*p2*a0 + 0.125;

        coef[3]=  0.125*p2*a1 + 0.5*CFL + 0.125;
        coef[4]=  -0.25*p2*a1 + 0.75;
        coef[5]=  0.125*p2*a1 - 0.5*CFL + 0.125;
        break;
    case 2: //LSP - central time - Peicewise parabolic
//
//    [-C^2*a^2 + 1/10, 2*C^2*a^2 + 4/5, -C^2*a^2 + 1/10, 0, 0]
        coef[0]= -a0 + 0.1;
        coef[1]= 2.0*a0 + 0.8;
        coef[2]= -a0 + 0.1;
        coef[3] = 0.0;
        coef[4] = 0.0;

//
//[-C^2*a*b - 1/2*C*a - 1/2*C*b - 1/10, 2*C^2*a*b - 4/5, -C^2*a*b +1/2*C*a + 1/2*C*b - 1/10, 0, 0]
//        negate
        coef[5] = a1 + 0.5*CFL + 0.1;
        coef[6] = -2.0*a1 + 0.8;
        coef[7] = a1 - 0.5*CFL + 0.1;
        coef[8] = 0.0;
        coef[9] = 0.0;

//
//    [1/4*C^2*a^2 - 1/10, -2*C^2*a^2 + 1/5, 7/2*C^2*a^2 + 4/5, -2*C^2*a^2 +1/5, 1/4*C^2*a^2 - 1/10]

        coef[10] = 0.25*a0 - 0.1;
        coef[11] = -2.0*a0 + 0.2;
        coef[12] = 3.5*a0 + 0.8;
        coef[13] = -2.0*a0 + 0.2;
        coef[14] = 0.25*a0 - 0.1;
//
//
//   [1/4*C^2*a*b + 1/4*C*a + 1/4*C*b + 1/10,
//	    -2*C^2*a*b - C*a - C*b - 1/5,
//	7/2*C^2*a*b - 4/5,
//	-2*C^2*a*b + C*a + C*b - 1/5,
//	1/4*C^2*a*b - 1/4*C*a -1/4*C*b + 1/10]
//	    negate
        coef[15] = -0.25*a1 - 0.25*CFL - 0.1;
        coef[16] = 2.0*a1 + CFL + 0.2;
        coef[17] = -3.5*a1 + 0.8;
        coef[18] = 2*a1 - CFL + 0.2;
        coef[19] = -0.25*a1 + 0.25*CFL - 0.1;

        break;
    case 3: //LSP - Fourier 4 - linear window
        double p3,den0,den1,den2,den3;
        p3=p2*pi;
        den0=(4.0/p2 + 1.0);
        den1=(30.0/p2 + 1.0);
        den2=(14.0/p2 + 3.0);
        den3=(32.0/p2 + 9.0);

        int ic=0;
        coef[ic++] = -(0.25*a0 - 1.0/6.0)/den0;
        coef[ic++] = -0.25*a0/den0;
        coef[ic++] = 4.0/(den0*p3);
        coef[ic++] = 0.0;
        coef[ic++] = (0.5*a0 + 2.0/3.0)/den0;
        coef[ic++] = (0.5*a0 + 4.0/p2)/den0;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = -(0.25*a0 - 1.0/6.0)/den0;
        coef[ic++] = -0.25*a0/den0;
        coef[ic++] = -4.0/(den0*p3);
        coef[ic++] = 0.0;

        coef[ic++] = (0.25*a1 + 0.25*CFL + 1.0/6.0)/den0;
        coef[ic++] = (0.25*a1 + CFL/p2)/den0;
        coef[ic++] = (0.5*CFL/pi + 4.0/p3)/den0;
        coef[ic++] = 0.25*(CFL/pi)/den0;
        coef[ic++] = -(0.5*a1 - 2.0/3.0)/den0;
        coef[ic++] = -(0.5*a1 - 4.0/p2)/den0;
        coef[ic++] = -(CFL/pi)/den0;
        coef[ic++] = -0.5*(CFL/pi)/den0;
        coef[ic++] = (0.25*a1 - 0.25*CFL + 1.0/6.0)/den0;
        coef[ic++] = (0.25*a1 - CFL/p2)/den0;
        coef[ic++] = (0.5*CFL/pi - 4.0/p3)/den0;
        coef[ic++] = 0.25*(CFL/pi)/den0;


        coef[ic++] = -1.5*a0/den1;
        coef[ic++] = -0.125*(p2*a0 - 9.0*a0 - 12.0/p2 + 4.0)/den1;
        coef[ic++] = -0.75*pi*a0/den1;
        coef[ic++] = -4.0/9.0*(13.0*a0/pi + 28.0/p3)/den1;
        coef[ic++] = 3.0*(a0 + 8.0/p2)/den1;
        coef[ic++] = 0.25*(2.0*p2*a0 + 9.0*a0 + 12.0/p2 + 8.0)/den1;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = -1.5*a0/den1;
        coef[ic++] = -0.125*(p2*a0 - 9.0*a0 - 12.0/p2 + 4.0)/den1;
        coef[ic++] = 0.75*pi*a0/den1;
        coef[ic++] = 4.0/9.0*(13.0*a0/pi + 28.0/p3)/den1;

        coef[ic++] = 1.5*(a1 + 4.0*CFL/p2)/den1;
        coef[ic++] = 0.125*(p2*a1 - 9.0*a1 - 6.0*CFL + 12.0/p2 - 4.0)/den1;
        coef[ic++] = 0.25*(3.0*pi*a1 + pi*CFL)/den1;
        coef[ic++] = 2.0/9.0*(26.0*a1/pi - 9.0*CFL/pi - 56.0/p3)/den1;
        coef[ic++] = -3.0*(a1 - 8.0/p2)/den1;
        coef[ic++] = -0.25*(2.0*p2*a1 + 9.0*a1 - 12.0/p2 - 8.0)/den1;
        coef[ic++] = -(pi*CFL)/den1;
        coef[ic++] = -28.0/3.0*(CFL/pi)/den1;
        coef[ic++] = 1.5*(a1 - 4.0*CFL/p2)/den1;
        coef[ic++] = 0.125*(p2*a1 - 9.0*a1 + 6.0*CFL + 12.0/p2 - 4.0)/den1;
        coef[ic++] = -0.25*(3.0*pi*a1 - pi*CFL)/den1;
        coef[ic++] = -2.0/9.0*(26.0*a1/pi + 9.0*CFL/pi - 56.0/p3)/den1;


        coef[ic++] = -72.0/(den2*p3);
        coef[ic++] = 9.0/4.0*pi*a0/den2;
        coef[ic++] = -0.375*(p2*a0 - 3.0*a0 + 12.0/p2 + 4.0)/den2;
        coef[ic++] = 6.0*a0/den2;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = 0.75*(2.0*p2*a0 + 3.0*a0 - 12.0/p2 + 8.0)/den2;
        coef[ic++] = 4.0*(5.0*a0 + 8.0/p2)/den2;
        coef[ic++] = 72.0/(den2*p3);
        coef[ic++] = -9.0/4.0*pi*a0/den2;
        coef[ic++] = -0.375*(p2*a0 - 3.0*a0 + 12.0/p2 + 4.0)/den2;
        coef[ic++] = 6.0*a0/den2;

        coef[ic++] = -9.0*(CFL/pi + 8.0/p3)/den2;
        coef[ic++] = -0.75*(3.0*pi*a1 + pi*CFL)/den2;
        coef[ic++] = 0.375*(p2*a1 - 3.0*a1 - 6.0*CFL - 12.0/p2 - 4.0)/den2;
        coef[ic++] = -2.0/3.0*(9.0*a1 + 40.0*CFL/p2)/den2;
        coef[ic++] = 18.0*(CFL/pi)/den2;
        coef[ic++] = 3.0*(pi*CFL)/den2;
        coef[ic++] = -0.75*(2.0*p2*a1 + 3.0*a1 + 12.0/p2 - 8.0)/den2;
        coef[ic++] = -4.0*(5.0*a1 - 8.0/p2)/den2;
        coef[ic++] = -9.0*(CFL/pi - 8.0/p3)/den2;
        coef[ic++] = 0.75*(3.0*pi*a1 - pi*CFL)/den2;
        coef[ic++] = 0.375*(p2*a1 - 3.0*a1 + 6.0*CFL - 12.0/p2 - 4.0)/den2;
        coef[ic++] = -2.0/3.0*(9.0*a1 - 40.0*CFL/p2)/den2;


        coef[ic++] = 0.0;
        coef[ic++] = 4.0/3.0*(13.0*a0/pi + 28.0/p3)/den3;
        coef[ic++] = 6.0*a0/den3;
        coef[ic++] = 0.375*(4.0*p2*a0 - 3.0*a0 + 3.0/p2 + 4.0)/den3;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = 4.0*(5.0*a0 + 8.0/p2)/den3;
        coef[ic++] = 0.75*(8.0*p2*a0 + 3.0*a0 - 3.0/p2 + 8.0)/den3;
        coef[ic++] = 0.0;
        coef[ic++] = -4.0/3.0*(13.0*a0/pi + 28.0/p3)/den3;
        coef[ic++] = 6.0*a0/den3;
        coef[ic++] = 0.375*(4.0*p2*a0 - 3.0*a0 + 3.0/p2 + 4.0)/den3;

        coef[ic++] = -4.5*(CFL/pi)/den3;
        coef[ic++] = -2.0/3.0*(26.0*a1/pi - 9.0*CFL/pi - 56.0/p3)/den3;
        coef[ic++] = -2.0/3.0*(9.0*a1 + 40.0*CFL/p2)/den3;
        coef[ic++] = -0.375*(4.0*p2*a1 - 3.0*a1 - 6.0*CFL - 3.0/p2 - 4.0)/den3;
        coef[ic++] = 9.0*(CFL/pi)/den3;
        coef[ic++] = 28.0*(CFL/pi)/den3;
        coef[ic++] = -4.0*(5.0*a1 - 8.0/p2)/den3;
        coef[ic++] = -0.75*(8.0*p2*a1 + 3.0*a1 + 3.0/p2 - 8.0)/den3;
        coef[ic++] = -4.5*(CFL/pi)/den3;
        coef[ic++] = 2.0/3.0*(26.0*a1/pi + 9.0*CFL/pi - 56.0/p3)/den3;
        coef[ic++] = -2.0/3.0*(9.0*a1 - 40.0*CFL/p2)/den3;
        coef[ic++] = -0.375*(4.0*p2*a1 - 3.0*a1 + 6.0*CFL - 3.0/p2 - 4.0)/den3;

        setupTrans();
    }
}


void ImpWidget::lspa() {
    double a0,a1,a2,p2,cib,cb;
    a0 = impl * impl * CFL * CFL;//a2
    a1 = impl * beta * CFL * CFL;//ab
    a2 = impl * back * CFL * CFL;//aB
    cib = CFL *( impl+beta );
    cb = CFL*back;
    p2 = pi * pi;

    switch(ibase) {
    default:
    case 0: // LSP - Adams - linear
        //[-C^2*a^2 + 1/6, 2*C^2*a^2 + 2/3, -C^2*a^2 + 1/6
        // -C^2*a*b - 1/2*C*a - 1/2*C*b - 1/6, 2*C^2*a*b - 2/3, -C^2*a*b + 1/2*C*a + 1/2*C*b - 1/6
        // -C^2*a*c - 1/2*C*c, 2*C^2*a*c, -C^2*a*c + 1/2*C*c]
        coef[0] = 1.0 / 6.0 - a0;
        coef[1] = 2.0 * ( 1.0 / 3.0 + a0 );
        coef[2] = 1.0 / 6.0 - a0;

        coef[3] = ( 1.0 / 3.0 + cib ) / 2.0 + a1;
        coef[4] = 2.0 * ( 1.0 / 3.0 - a1 );
        coef[5] = ( 1.0 / 3.0 - cib ) / 2.0 + a1;

        coef[6] = cb/2.0 + a2;
        coef[7] = -2.0*a2;
        coef[8] = -cb/2.0 + a2;

        break;
    case 1: // LSP - Adams -raised cosine
        //[-1/8*pi^2*C^2*a^2 + 1/8,   1/4*pi^2*C^2*a^2 + 3/4,   -1/8*pi^2*C^2*a^2 +1/8,
        //-1/8*pi^2*C^2*a*b - 1/2*C*a - 1/2*C*b - 1/8,   1/4*pi^2*C^2*a*b -3/4,   -1/8*pi^2*C^2*a*b + 1/2*C*a + 1/2*C*b - 1/8,
        //  -1/8*pi^2*C^2*a*c -1/2*C*c, 1/4*pi^2*C^2*a*c, -1/8*pi^2*C^2*a*c + 1/2*C*c

        coef[0] = ( 1.0 - a0 * p2 ) / 8.0;
        coef[1] = ( 3.0 + a0 * p2 ) / 4.0;
        coef[2] = ( 1.0 - a0 * p2 ) / 8.0;

        coef[3] = ( 1.0 + 4.0 * cib + a1 * p2 ) / 8.0;
        coef[4] = ( 3.0 - a1 * p2 ) / 4.0;
        coef[5] = ( 1.0 - 4.0 * cib + a1 * p2 ) / 8.0;

        coef[6] = cb/2.0 + a2*p2/8.0;
        coef[7] = -a2*p2/4.0;
        coef[8] = -cb/2.0 + a2*p2/8.0;
        break;
    case 2: // LSP - Adams - Peicewise parabolic
        /*[-C^2*a^2 + 1/10, 2*C^2*a^2 + 4/5, -C^2*a^2 + 1/10, 0, 0]
        */
        coef[0] = -a0 + 0.1;
        coef[1] =  2.0*a0 + 0.8;
        coef[2] =  -a0 + 0.1;
        coef[3] =  0;
        coef[4] =  0;

//  [-C^2*a*b - 1/2*C*a - 1/2*C*b - 1/10, 2*C^2*a*b - 4/5, -C^2*a*b + 1/2*C*a + 1/2*C*b - 1/10, 0, 0]
//       NEGATE
        coef[5] = a1 + 0.5*cib + 0.1;
        coef[6] = -2.0*a1 + 0.8;
        coef[7] = a1 - 0.5*cib + 0.1;
        coef[8] =  0;
        coef[9] =  0;

//
//	[-C^2*a*c - 1/2*C*c, 2*C^2*a*c, -C^2*a*c + 1/2*C*c, 0, 0]
//	    NEGATE
        coef[10] = -a2 - 0.5*cb;
        coef[11] =  2.0*a2;
        coef[12] = -a2 + 0.5*cb;
        coef[13] =  0;
        coef[14] =  0;

//
//	 [1/4*C^2*a^2 - 1/10, -2*C^2*a^2 + 1/5, 7/2*C^2*a^2 + 4/5, -2*C^2*a^2 +1/5, 1/4*C^2*a^2 - 1/10]
        coef[15] = 0.25*a0 - 0.1;
        coef[16] = -2.0*a0 + 0.2;
        coef[17] = 3.5*a0 + 0.8;
        coef[18] = -2.0*a0 + 0.2;
        coef[19] = 0.25*a0 - 0.1;

//
//	[1/4*C^2*a*b + 1/4*C*a + 1/4*C*b + 1/10,
//	    -2*C^2*a*b - C*a - C*b - 1/5,
//	7/2*C^2*a*b - 4/5,
//	    -2*C^2*a*b + C*a + C*b - 1/5,
//	    1/4*C^2*a*b - 1/4*C*a - 1/4*C*b + 1/10]
//	    NEGATE
        coef[20] =  -0.25*a1 - 0.25*cib - 0.1;
        coef[21] =  2.0*a1 + cib + 0.2;
        coef[22] =  -3.5*a1 + 0.8;
        coef[23] =  2.0*a1 - cib + 0.2;
        coef[24] =  -0.25*a1 + 0.25*cib - 0.1;

//	    [1/4*C^2*a*c + 1/4*C*c,
//	    -2*C^2*a*c - C*c,
//	    7/2*C^2*a*c,
//	    -2*C^2*a*c + C*c,
//	1/4*C^2*a*c - 1/4*C*c]
//	     NEGATE
        coef[25] =  -0.25*a2 - 0.25*cb;
        coef[26] =  2.0*a2 + cb;
        coef[27] =  -3.5*a2;
        coef[28] =  2.0*a2 - cb;
        coef[29] =  -0.25*a2 + 0.25*cb;

        break;
    }
}


void ImpWidget::femc() {
    switch(ibase) {
    default:
        ibase = 0;
    case 0: // FEM - central time - linear
        coef[0] = -0.5*CFL*impl + 1.0/6.0;
        coef[1] =  2.0/3.0;
        coef[2] =  0.5*CFL*impl + 1.0/6.0;

        coef[3] =  0.5*CFL*beta + 1.0/6.0;
        coef[4] =  2.0/3.0;
        coef[5] =  -0.5*CFL*beta + 1.0/6.0;
        break;
    case 1: // FEM - central time - raised cosine
        /*[-1/2*C*a + 1/8, 3/4, 1/2*C*a + 1/8, -1/2*C*b - 1/8, -3/4, 1/2*C*b -1/8]*/
        coef[0] = -0.5*CFL*impl + 0.125;
        coef[1] =  0.75;
        coef[2] =  0.5*CFL*impl + 0.125;

        coef[3] =  0.5*CFL*beta + 0.125;
        coef[4] =  0.75;
        coef[5] =  -0.5*CFL*beta + 0.125;
        break;
    case 2: // FEM - central time - Peicewise parabolic
//
//   [-1/2*C*a + 1/10, 4/5, 1/2*C*a + 1/10, 0, 0]
        coef[0] = -0.5*CFL*impl + 0.1;
        coef[1] =  0.8;
        coef[2] =  0.5*CFL*impl + 0.1;
        coef[3] =  0.0;
        coef[4] =  0.0;
//
//   [-1/2*C*b - 1/10, -4/5, 1/2*C*b - 1/10, 0, 0]
// negate
        coef[5] =  0.5*CFL*beta + 0.1;
        coef[6] =  0.8;
        coef[7] =  -0.5*CFL*beta + 0.1;
        coef[8] =  0.0;
        coef[9] =  0.0;
//
//    [1/4*C*a - 1/10, -C*a + 1/5, 4/5, C*a + 1/5, -1/4*C*a - 1/10]
        coef[10] = 0.25*CFL*impl - 0.1;
        coef[11] = -CFL*impl + 0.2;
        coef[12] = 0.8;
        coef[13] = CFL*impl + 0.2;
        coef[14] = -0.25*CFL*impl - 0.1;
//
//    [1/4*C*b + 1/10, -C*b - 1/5, -4/5, C*b - 1/5, -1/4*C*b + 1/10]
// negate

        coef[15] = -0.25*CFL*beta - 0.1;
        coef[16] = CFL*beta + 0.2;
        coef[17] = .8;
        coef[18] = -CFL*beta + 0.2;
        coef[19] = 0.25*CFL*beta - 0.1;
        break;
    case 3: // FEM - central time - Fourier 4 Linear
        double p2,p3,p4,den0,denf1,denf2,denf3,cimp,cbet;
	int ic;
	cimp = CFL*impl;
	cbet = CFL*beta;
	p2=pi*pi;
        p3=p2*pi;
	p4=p3*pi;
        den0=(4.0/p2 + 1.0);
        denf1=((3*pi - 2*p3)/p3 + 57/p2 + 4);
        denf2=(3*(3*pi + 2*p3)/p3 - 37/p2 - 12);
        denf3=(3*(3*pi + 8*p3)/p3 + 247/p2 + 48);
        ic=0;

        coef[ic++] = -1.0/12.0*(3*p3*cimp - 2*p3)/(den0*p3);
        coef[ic++] = 0.5*((2*pi - pi*cimp)/p3 - (2*pi + pi*cimp)/p3)/den0;
        coef[ic++] = -0.5*((p2*cimp - 4)/p3 - 4/p3)/den0;
        coef[ic++] = -0.25*((p2*cimp - 1)/p3 + 1/p3)/den0;
        coef[ic++] = -1.0/12.0*((3*p3*cimp - 4*p3)/p3 - (3*p3*cimp + 4*p3)/p3)/den0;
        coef[ic++] = 0.5*((4*pi - (pi - p3)*cimp)/p3 + (4*pi + (pi - p3)*cimp)/p3)/den0;
        coef[ic++] = 0.5*((p2*cimp - 2*p2 + 4)/p3 + (p2*cimp + 2*p2 - 4)/p3)/den0;
        coef[ic++] = 0.25*((p2*cimp - 2*p2 + 1)/p3 + (p2*cimp + 2*p2 - 1)/p3)/den0;
        coef[ic++] = 1.0/12.0*(3*p3*cimp + 2*p3)/(den0*p3);
        coef[ic++] = -0.5*((2*pi - pi*cimp)/p3 - (2*pi + pi*cimp)/p3)/den0;
        coef[ic++] = -0.5*((p2*cimp + 4)/p3 + 4/p3)/den0;
        coef[ic++] = -0.25*((p2*cimp + 1)/p3 - 1/p3)/den0;

        coef[ic++] = 1.0/12.0*(3*p3*cbet + 2*p3)/(den0*p3);
        coef[ic++] = -0.5*((2*pi - pi*cbet)/p3 - (2*pi + pi*cbet)/p3)/den0;
        coef[ic++] = 0.5*((p2*cbet + 4)/p3 + 4/p3)/den0;
        coef[ic++] = 0.25*((p2*cbet + 1)/p3 - 1/p3)/den0;
        coef[ic++] = -1.0/12.0*((3*p3*cbet - 4*p3)/p3 - (3*p3*cbet + 4*p3)/p3)/den0;
        coef[ic++] = 0.5*((4*pi - (pi - p3)*cbet)/p3 + (4*pi + (pi - p3)*cbet)/p3)/den0;
        coef[ic++] = -0.5*((p2*cbet - 2*p2 + 4)/p3 + (p2*cbet + 2*p2 - 4)/p3)/den0;
        coef[ic++] = -0.25*((p2*cbet - 2*p2 + 1)/p3 + (p2*cbet + 2*p2 - 1)/p3)/den0;
        coef[ic++] = -1.0/12.0*(3*p3*cbet - 2*p3)/(den0*p3);
        coef[ic++] = 0.5*((2*pi - pi*cbet)/p3 - (2*pi + pi*cbet)/p3)/den0;
        coef[ic++] = 0.5*((p2*cbet - 4)/p3 - 4/p3)/den0;
        coef[ic++] = 0.25*((p2*cbet - 1)/p3 + 1/p3)/den0;


        coef[ic++] = 6*((2*pi - pi*cimp)/p3 - (2*pi + pi*cimp)/p3)/denf1;
        coef[ic++] = 0.5*((3*pi + 3*p3*cimp - 2*p3)/p3 + 3/p2)/denf1;
        coef[ic++] = -0.25*(((2*p4 - 9*p2)*cimp + 6)/p3 + 3*(3*p2*cimp - 2)/p3)/denf1;
        coef[ic++] = -4.0/9.0*(2*(3*p2*cimp + 14)/p3 - (15*p2*cimp - 28)/p3)/denf1;
        coef[ic++] = 6*((4*pi - pi*cimp)/p3 + (4*pi + pi*cimp)/p3)/denf1;
        coef[ic++] = (3*(pi - p3*cimp)/p3 + 3*(pi + p3*cimp)/p3 + 4)/denf1;
        coef[ic++] = 0.5*(((2*p4 + 3*p2)*cimp - 6*p2 + 3)/p3 + ((2*p4 + 3*p2)*cimp + 6*p2 - 3)/p3 - 3*(p2*cimp - 4*p2 + 1)/p3 - 3*(p2*cimp + 4*p2 - 1)/p3)/denf1;
        coef[ic++] = 4.0/9.0*((21*p2*cimp - 18*p2 + 28)/p3 + (21*p2*cimp + 18*p2 - 28)/p3)/denf1;
        coef[ic++] = -6*((2*pi - pi*cimp)/p3 - (2*pi + pi*cimp)/p3)/denf1;
        coef[ic++] = 0.5*((3*pi - 3*p3*cimp - 2*p3)/p3 + 3/p2)/denf1;
        coef[ic++] = -0.25*(((2*p4 - 9*p2)*cimp - 6)/p3 + 3*(3*p2*cimp + 2)/p3)/denf1;
        coef[ic++] = -4.0/9.0*(2*(3*p2*cimp - 14)/p3 - (15*p2*cimp + 28)/p3)/denf1;

        coef[ic++] = -6*((2*pi - pi*cbet)/p3 - (2*pi + pi*cbet)/p3)/denf1;
        coef[ic++] = 0.5*((3*pi - 3*p3*cbet - 2*p3)/p3 + 3/p2)/denf1;
        coef[ic++] = 0.25*(((2*p4 - 9*p2)*cbet - 6)/p3 + 3*(3*p2*cbet + 2)/p3)/denf1;
        coef[ic++] = 4.0/9.0*(2*(3*p2*cbet - 14)/p3 - (15*p2*cbet + 28)/p3)/denf1;
        coef[ic++] = 6*((4*pi - pi*cbet)/p3 + (4*pi + pi*cbet)/p3)/denf1;
        coef[ic++] = (3*(pi - p3*cbet)/p3 + 3*(pi + p3*cbet)/p3 + 4)/denf1;
        coef[ic++] = -0.5*(((2*p4 + 3*p2)*cbet - 6*p2 + 3)/p3 + ((2*p4 + 3*p2)*cbet + 6*p2 - 3)/p3 - 3*(p2*cbet - 4*p2 + 1)/p3 - 3*(p2*cbet + 4*p2 - 1)/p3)/denf1;
        coef[ic++] = -4.0/9.0*((21*p2*cbet - 18*p2 + 28)/p3 + (21*p2*cbet + 18*p2 - 28)/p3)/denf1;
        coef[ic++] = 6*((2*pi - pi*cbet)/p3 - (2*pi + pi*cbet)/p3)/denf1;
        coef[ic++] = 0.5*((3*pi + 3*p3*cbet - 2*p3)/p3 + 3/p2)/denf1;
        coef[ic++] = 0.25*(((2*p4 - 9*p2)*cbet + 6)/p3 + 3*(3*p2*cbet - 2)/p3)/denf1;
        coef[ic++] = 4.0/9.0*(2*(3*p2*cbet + 14)/p3 - (15*p2*cbet - 28)/p3)/denf1;


        coef[ic++] = -18*((p2*cimp - 4)/p3 - 4/p3)/denf2;
        coef[ic++] = -0.75*(((2*p4 + 9*p2)*cimp - 6)/p3 - 3*(3*p2*cimp - 2)/p3)/denf2;
        coef[ic++] = 1.5*((3*pi - 3*p3*cimp + 2*p3)/p3 + 3/p2)/denf2;
        coef[ic++] = 16.0/3.0*((3*pi - 5*pi*cimp)/p3 - (3*pi + 5*pi*cimp)/p3)/denf2;
        coef[ic++] = 18*((p2*cimp - 2*p2 + 4)/p3 + (p2*cimp + 2*p2 - 4)/p3)/denf2;
        coef[ic++] = 1.5*(((2*p4 - 3*p2)*cimp - 6*p2 + 3)/p3 + ((2*p4 - 3*p2)*cimp + 6*p2 - 3)/p3 + 3*(p2*cimp - 4*p2 + 1)/p3 + 3*(p2*cimp + 4*p2 - 1)/p3)/denf2;
        coef[ic++] = 3*((3*p3*cimp - 2*p3)/p3 - (3*p3*cimp + 2*p3)/p3 + 3*(pi - p3*cimp)/p3 + 3*(pi + p3*cimp)/p3)/denf2;
        coef[ic++] = -4.0/3.0*((24*pi - (20*pi - 9*p3)*cimp)/p3 + (24*pi + (20*pi - 9*p3)*cimp)/p3)/denf2;
        coef[ic++] = -18*((p2*cimp + 4)/p3 + 4/p3)/denf2;
        coef[ic++] = -0.75*(((2*p4 + 9*p2)*cimp + 6)/p3 - 3*(3*p2*cimp + 2)/p3)/denf2;
        coef[ic++] = 1.5*((3*pi + 3*p3*cimp + 2*p3)/p3 + 3/p2)/denf2;
        coef[ic++] = -16.0/3.0*((3*pi - 5*pi*cimp)/p3 - (3*pi + 5*pi*cimp)/p3)/denf2;

        coef[ic++] = 18*((p2*cbet + 4)/p3 + 4/p3)/denf2;
        coef[ic++] = 0.75*(((2*p4 + 9*p2)*cbet + 6)/p3 - 3*(3*p2*cbet + 2)/p3)/denf2;
        coef[ic++] = 1.5*((3*pi + 3*p3*cbet + 2*p3)/p3 + 3/p2)/denf2;
        coef[ic++] = -16.0/3.0*((3*pi - 5*pi*cbet)/p3 - (3*pi + 5*pi*cbet)/p3)/denf2;
        coef[ic++] = -18*((p2*cbet - 2*p2 + 4)/p3 + (p2*cbet + 2*p2 - 4)/p3)/denf2;
        coef[ic++] = -1.5*(((2*p4 - 3*p2)*cbet - 6*p2 + 3)/p3 + ((2*p4 - 3*p2)*cbet + 6*p2 - 3)/p3 + 3*(p2*cbet - 4*p2 + 1)/p3 + 3*(p2*cbet + 4*p2 - 1)/p3)/denf2;
        coef[ic++] = 3*((3*p3*cbet - 2*p3)/p3 - (3*p3*cbet + 2*p3)/p3 + 3*(pi - p3*cbet)/p3 + 3*(pi + p3*cbet)/p3)/denf2;
        coef[ic++] = -4.0/3.0*((24*pi - (20*pi - 9*p3)*cbet)/p3 + (24*pi + (20*pi - 9*p3)*cbet)/p3)/denf2;
        coef[ic++] = 18*((p2*cbet - 4)/p3 - 4/p3)/denf2;
        coef[ic++] = 0.75*(((2*p4 + 9*p2)*cbet - 6)/p3 - 3*(3*p2*cbet - 2)/p3)/denf2;
        coef[ic++] = 1.5*((3*pi - 3*p3*cbet + 2*p3)/p3 + 3/p2)/denf2;
        coef[ic++] = 16.0/3.0*((3*pi - 5*pi*cbet)/p3 - (3*pi + 5*pi*cbet)/p3)/denf2;


        coef[ic++] = 36*((p2*cimp - 1)/p3 + 1/p3)/denf3;
        coef[ic++] = 16.0/3.0*(2*(3*p2*cimp + 14)/p3 - (15*p2*cimp - 28)/p3)/denf3;
        coef[ic++] = -64.0/3.0*((3*pi - 5*pi*cimp)/p3 - (3*pi + 5*pi*cimp)/p3)/denf3;
        coef[ic++] = 1.5*((3*pi - 12*p3*cimp + 8*p3)/p3 + 3/p2)/denf3;
        coef[ic++] = -36*((p2*cimp - 2*p2 + 1)/p3 + (p2*cimp + 2*p2 - 1)/p3)/denf3;
        coef[ic++] = -16.0/3.0*((21*p2*cimp - 18*p2 + 28)/p3 + (21*p2*cimp + 18*p2 - 28)/p3)/denf3;
        coef[ic++] = 16.0/3.0*((24*pi - (20*pi - 9*p3)*cimp)/p3 + (24*pi + (20*pi - 9*p3)*cimp)/p3)/denf3;
        coef[ic++] = -3*(4*(3*p3*cimp - 2*p3)/p3 - 4*(3*p3*cimp + 2*p3)/p3 + 3*(pi - 4*p3*cimp)/p3 + 3*(pi + 4*p3*cimp)/p3)/denf3;
        coef[ic++] = 36*((p2*cimp + 1)/p3 - 1/p3)/denf3;
        coef[ic++] = 16.0/3.0*(2*(3*p2*cimp - 14)/p3 - (15*p2*cimp + 28)/p3)/denf3;
        coef[ic++] = 64.0/3.0*((3*pi - 5*pi*cimp)/p3 - (3*pi + 5*pi*cimp)/p3)/denf3;
        coef[ic++] = 1.5*((3*pi + 12*p3*cimp + 8*p3)/p3 + 3/p2)/denf3;

        coef[ic++] = -36*((p2*cbet + 1)/p3 - 1/p3)/denf3;
        coef[ic++] = -16.0/3.0*(2*(3*p2*cbet - 14)/p3 - (15*p2*cbet + 28)/p3)/denf3;
        coef[ic++] = 64.0/3.0*((3*pi - 5*pi*cbet)/p3 - (3*pi + 5*pi*cbet)/p3)/denf3;
        coef[ic++] = 1.5*((3*pi + 12*p3*cbet + 8*p3)/p3 + 3/p2)/denf3;
        coef[ic++] = 36*((p2*cbet - 2*p2 + 1)/p3 + (p2*cbet + 2*p2 - 1)/p3)/denf3;
        coef[ic++] = 16.0/3.0*((21*p2*cbet - 18*p2 + 28)/p3 + (21*p2*cbet + 18*p2 - 28)/p3)/denf3;
        coef[ic++] = 16.0/3.0*((24*pi - (20*pi - 9*p3)*cbet)/p3 + (24*pi + (20*pi - 9*p3)*cbet)/p3)/denf3;
        coef[ic++] = -3*(4*(3*p3*cbet - 2*p3)/p3 - 4*(3*p3*cbet + 2*p3)/p3 + 3*(pi - 4*p3*cbet)/p3 + 3*(pi + 4*p3*cbet)/p3)/denf3;
        coef[ic++] = -36*((p2*cbet - 1)/p3 + 1/p3)/denf3;
        coef[ic++] = -16.0/3.0*(2*(3*p2*cbet + 14)/p3 - (15*p2*cbet - 28)/p3)/denf3;
        coef[ic++] = -64.0/3.0*((3*pi - 5*pi*cbet)/p3 - (3*pi + 5*pi*cbet)/p3)/denf3;
        coef[ic++] = 1.5*((3*pi - 12*p3*cbet + 8*p3)/p3 + 3/p2)/denf3;
    }
}

void ImpWidget::fema() {
    switch(ibase) {
    default:
    case 0: // linear
    case 1: // raised cosine
        ibase = 2;
    case 2: // Peicewise parabolic
//
//[-1/2*C*a + 1/10, 4/5, 1/2*C*a + 1/10, 0, 0]
        coef[0]= -0.5*CFL*impl +0.1;
        coef[1]= 0.8;
        coef[2]= 0.5*CFL*impl + 0.1;
        coef[3]= 0;
        coef[4]= 0;
//
//[-1/2*C*b - 1/10, -4/5, 1/2*C*b - 1/10, 0, 0]
// NEGATE
        coef[5]= 0.5*CFL*beta + 0.1;
        coef[6]= 0.8;
        coef[7]= -0.5*CFL*beta + 0.1;
        coef[8]= 0;
        coef[9]= 0;
//
//[-1/2*C*c, 0, 1/2*C*c, 0, 0]
//	  0]  NEGATE
        coef[10]= 0.5*CFL*back;
        coef[11]= 0.0;
        coef[12]= -0.5*CFL*back;
        coef[13]= 0;
        coef[14]= 0;
//
//
// [1/4*C*a - 1/10, -C*a + 1/5, 4/5, C*a + 1/5, -1/4*C*a - 1/10]
        coef[15]= 0.25*CFL*impl - 0.1 ;
        coef[16]= -CFL*impl + 0.2;
        coef[17]= 0.8;
        coef[18]= CFL*impl + 0.2;
        coef[19]= -0.25*CFL*impl -0.1;
//
//[1/4*C*b + 1/10, -C*b - 1/5, -4/5, C*b - 1/5, -1/4*C*b + 1/10]
//	      NEGATE
        coef[20]= -0.25*CFL*beta - 0.1 ;
        coef[21]= CFL*beta + 0.2;
        coef[22]= 0.8;
        coef[23]= -CFL*beta + 0.2;
        coef[24]= 0.25*CFL*beta - 0.1;
//
//[1/4*C*c, -C*c, 0, C*c, -1/4*C*c]
//  NEGATE
        coef[25]= -0.25*CFL*back;
        coef[26]= CFL*back;
        coef[27]= 0.0;
        coef[28]= -CFL*back;
        coef[29]= 0.25*CFL*back;
//
//
        break;
    }
}

void ImpWidget::initSin(const double value) {
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
    if(transform) forwardTrans(U_,Utran);
    unstable = false;
}

double* ImpWidget::getU() {
    if(transform) {
        reversTrans(Utran,U_);
    }
    return U_;
}


void ImpWidget::reversTrans(double *from,double *to) {
    int nn,nb;
    for(nn=0; nn<N_ ; nn += nblock) {
        for( nb=0; nb < nblock; nb++) {
            gsl_vector_set(TranVec,nb,from[nn+nb]);
        }
        gsl_blas_dgemv(CblasNoTrans,1.0,Mback,TranVec,0.0,UVec);
        for( nb=0; nb < nblock; nb++) {
            to[nn+nb] = gsl_vector_get(UVec,nb);
        }
    }
}

void ImpWidget::forwardTrans(double *from,double *to) {
    int nn,nb;
    for(nn=0; nn<N_ ; nn += nblock) {
        for( nb=0; nb < nblock; nb++) {
            gsl_vector_set(UVec,nb,from[nn+nb]);
        }
        gsl_linalg_LU_solve(Mforw,permut,UVec,TranVec);
        for( nb=0; nb < nblock; nb++) {
            to[nn+nb] = gsl_vector_get(TranVec,nb);
        }
    }
}

void ImpWidget::setupTrans() {
    double testin[] = { 1.0, 2.0, 3.0, 4.0 };
    double testout[4];
    if(ibase == 3) {
        if( matSize_gsl != 4 ) {
            if( matSize_gsl ) {
                // deallocate matricies
                gsl_vector_free(UVec);
                gsl_vector_free(TranVec);
                gsl_matrix_free(Mforw);
                gsl_matrix_free(Mback);
            }
            UVec = gsl_vector_alloc(4);
            TranVec = gsl_vector_alloc(4);
            Mforw = gsl_matrix_alloc(4,4);
            Mback =  gsl_matrix_alloc(4,4);
            permut = gsl_permutation_alloc(4);
            matSize_gsl = 4;
        }
//At location  -3/2
//  f[4] =  1/4, f[5] =  -1/8*sqrt(2), f[6] =  -1/8*sqrt(2), f[7] =  1/4
        gsl_matrix_set (Mback,0,0,0.25);
        gsl_matrix_set (Mback,0,1,-sqrt(2)/8);
        gsl_matrix_set (Mback,0,2,-sqrt(2)/8);
        gsl_matrix_set (Mback,0,3,0.25);
//At location  -1/2
//  f[4] =  3/4, f[5] =  3/8*sqrt(2),  f[6] =  -3/8*sqrt(2), f[7] =  -3/4
        gsl_matrix_set (Mback,1,0,0.75);
        gsl_matrix_set (Mback,1,1,sqrt(2)*0.375);
        gsl_matrix_set (Mback,1,2,-sqrt(2)*0.375);
        gsl_matrix_set (Mback,1,3,-0.75);
//At location  1/2
//  f[4] =  3/4, f[5] =  3/8*sqrt(2),  f[6] =  3/8*sqrt(2),  f[7] =  3/4
        gsl_matrix_set (Mback,2,0,0.75);
        gsl_matrix_set (Mback,2,1,sqrt(2)*0.375);
        gsl_matrix_set (Mback,2,2,sqrt(2)*0.375);
        gsl_matrix_set (Mback,2,3,0.75);
//At location  3/2
//  f[4] =  1/4, f[5] =  -1/8*sqrt(2), f[6] =  1/8*sqrt(2),  f[7] =  -1/4
        gsl_matrix_set (Mback,3,0,0.25);
        gsl_matrix_set (Mback,3,1,-sqrt(2)/8);
        gsl_matrix_set (Mback,3,2,sqrt(2)/8);
        gsl_matrix_set (Mback,3,3,-0.25);

        gsl_matrix_memcpy (Mforw,Mback);
        gsl_linalg_LU_decomp (Mforw, permut, &signum);
        transform = true;

        return;
        // debug information printout
        std::cout << " Mback\n";
        gsl_matrix_fprintf (stdout,Mback,"%g");
        std::cout << "\n\n Mforw\n";
        gsl_matrix_fprintf (stdout,Mforw,"%g");
        gsl_vector_view bvec = gsl_vector_view_array (testin, 4);
        gsl_vector_view xvec = gsl_vector_view_array (testout, 4);
        std::cout << "\n\n B\n";
        gsl_vector_fprintf (stdout,&(bvec.vector),"%g");
        gsl_linalg_LU_solve (Mforw, permut, &(bvec.vector), &(xvec.vector));
        std::cout << "\n\n Btran\n";
        gsl_vector_fprintf (stdout,&(xvec.vector),"%g");
        gsl_blas_dgemv(CblasNoTrans,1.0,Mback,&(xvec.vector),0.0,&(bvec.vector));
        std::cout << "\n\n B\n";
        gsl_vector_fprintf (stdout,&(bvec.vector),"%g");


    }
}
bool ImpWidget::canSolve(int equ) {
    return (equ == 0);
}
