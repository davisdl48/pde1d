/*
    Computes Runge-Kutta advection solution using Butcher Tablue
    and Finite Element / LS derivatives
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

#include "rkwidget.h"
#include <iostream>

RKWidget::RKWidget ( QWidget* parent ) : SolvWidget(parent)
{   int iwid = 4;
    lwork = 0;
    nrhs  = 1;
    nStage = 0;
    setTitle(tr("Runge Kutta"));

    methodLabel = new QLabel ( tr ( "Method" ) );
    methodBox = new QComboBox ( this );
    methodBox->addItem( tr("Mid-Point"));
    methodBox->addItem( tr("Trapezoid"));
    methodBox->addItem( tr("RK4"));
    methodBox->addItem( tr("Arb-RK-4"));
    connect ( methodBox, SIGNAL( activated(int) ), this, SLOT( setMethod(int) ) );
    verticalLayout->insertWidget ( iwid++, methodLabel );
    verticalLayout->insertWidget ( iwid++, methodBox );
    weightLabel = new QLabel ( tr ( "Basis Functions" ) );
    weightBox = new QComboBox ( this );
    weightBox->addItem ( tr ( "Linear Basis" ) );
    weightBox->addItem ( tr ( "Cosine Basis" ) );
    weightBox->addItem ( tr ( "Finite Diff" ) );
    weightBox->addItem ( tr ( "Finite Diff - 4th ord" ) );
    //weightBox->addItem ( tr ( "Parabolic Basis" ) );
    //weightBox->addItem ( tr ( "Fourier 4 Basis" ) );
    verticalLayout->insertWidget ( iwid++, weightLabel );
    verticalLayout->insertWidget ( iwid++, weightBox );
    connect ( weightBox, SIGNAL( activated(int) ), this, SLOT( setBasis(int) ) );
    plotNameEdit->setText ( title );

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
    setColor ( Qt::magenta );
    ncoef = 0;
    aexist = false;
    dirty = true;
    nblock = 1;
    setSize(100);
    setBasis(0);
    method=-1;
    setMethod(0);
  unstable = false;
    //set_default_options(&options);
}

RKWidget::RKWidget ( const RKWidget& other )
{

}

RKWidget::~RKWidget()
{
    if ( N_ != 0 ) {

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
}

RKWidget& RKWidget::operator= ( const RKWidget& other )
{
    return *this;
}

bool RKWidget::operator== ( const RKWidget& other ) const
{
    return (this == &other);
}


void RKWidget::setBasis ( int index )
{
    dirty = true;
    switch(index) {
    default:
        ibase = 0;
    case 0: // linear
    case 1: // raised cosine
    case 2: // Finite Difference
        nvar = 3;
        nblock = 1;
        ivar = 1;
        break;
    case 3:// Finite Difference 4th
        nvar = 5;
        nblock = 1;
        ivar=2;
	break;
    case 4: // peicewise parabolic
        nvar = 5;
        nblock = 2;
        ivar=1;
	break;
    case 5: // Fourier 4
        nvar = 12;
        nblock = 4;
        ivar = 4;
        /// todo -  setupTrans() and friends  -- error wrong approach -- new envwidget
        break;
    }
    ibase = index;
    ndn = ivar;
    nup = nvar - ivar;
    weightBox->setCurrentIndex ( ibase );
    updateCoef(method);
}

void RKWidget::setSize ( const size_t value )
{
    int add = value%nblock;
    size_t newval = value + add;
    if ( newval == N_ ) return;
    dirty = true;
    if ( N_ != 0 ) {
        if(nStage != 0) delete[] b_k;

        SUPERLU_FREE (rhsb);
        SUPERLU_FREE (rhsx);
        SUPERLU_FREE (etree);
        SUPERLU_FREE (perm_r);
        SUPERLU_FREE (perm_c);
        SUPERLU_FREE (R);
        SUPERLU_FREE (C);
        SUPERLU_FREE (ferr);
        SUPERLU_FREE (berr);

        if(aexist) {
            // ??? Destroy_CompCol_Matrix(&A);
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
    resize(value);
    if(nStage != 0) b_k = new double[N_*nStage];

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

void RKWidget::step ( const size_t nStep )
{
    DNformat *Bstore;
    DNformat *Xstore;
    double temp;
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
        for( int ns = 0; ns < nStage; ns++ ) {
            ///set B matrix
            Bstore= (DNformat *)B.Store;
            rhsb=(double*)Bstore->nzval;
            if (nup+ndn == 0) {
                // solve seperate independent blocks???
                std::cout << "Discontinuous Galerkin Not Implimented\n";
                return;
            } else {
                fillB(ns);
            }

            ///solve factored system
            StatInit(&stat);
            options.Fact = FACTORED;
            dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
                   &L, &Up, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
                   &mem_usage, &stat, &info);

            //if( n == 1 ) printf("Triangular solve: dgssvx() returns info %d\n", info);

            ///update U_
            if ( info == 0 || info == N_+1 ) {

                /* This is how you could access the solution matrix. */
                Xstore = (DNformat *)X.Store;
                rhsx = (double*)(Xstore->nzval);
                for ( size_t i = 0; i <  N_ ; i++ ) {
                    if(rhsx[i]> 1e16) unstable = true;
                    b_k[i+N_*ns]=rhsx[i];
                }
            } else {
                std::cout << "ERROR: Matrix Solution Failed   info = " << info << std::endl;
            }
            StatFree(&stat);
        }
        for(int j=0 ; j<nStage; j++) {
            temp=b_b[j];
            if( temp != 0 ) {
                for(size_t i = 0; i <  N_ ; i++ ) {
                    U_[i] += temp*b_k[i+j*N_];
                }
            }
        }
        cStep++;
        totCFL += CFL;
    }
}


void RKWidget::updateCoef(int value) {
    if( value != method ) setMethod(value);
    double p2;
    int ic;
    p2 = pi*pi;
    if(ncoef != nblock*nvar*3) {
        if( ncoef != 0 ) {
            delete[] coef;
        }
        ncoef = nblock*nvar*3;
        coef = new double[ncoef];
    }
    switch(ibase) {
    default:
        ibase=0;
    case 0: // LS-FEM - linear
        /*[1/6, 2/3, 1/6, 1/2*C, 0, -1/2*C]*/
        coef[0]= 1.0/6;
        coef[1]= 2.0/3;
        coef[2]= 1.0/6;

        coef[3]=  1.0/2*CFL*e_;
        coef[4]=  0;
        coef[5]=  -1.0/2*CFL*e_;

        coef[6]= CFL*visc_*d_/dx;
        coef[7]= -2*CFL*visc_*d_/dx;
        coef[8]= CFL*visc_*d_/dx;
        break;
    case 1: // LS-FEM - raised cosine
        /*[1/8, 3/4, 1/8, 1/2*C, 0, -1/2*C]*/

        coef[0]= 1.0/8;
        coef[1]= 3.0/4;
        coef[2]= 1.0/8;

        coef[3]=  1.0/2*CFL*e_;
        coef[4]=  0;
        coef[5]=  -1.0/2*CFL*e_;

        coef[6]= CFL*visc_*d_*p2/dx/8;
        coef[7]= -CFL*visc_*d_*p2/dx/4;
        coef[8]= CFL*visc_*d_*p2/dx/8;

        break;
    case 2: // Finite Difference/Volume
        /*[1/8, 3/4, 1/8, 1/2*C, 0, -1/2*C]*/

        coef[0]= 0;
        coef[1]= 1;
        coef[2]= 0;

        coef[3]=  1.0/2*CFL*e_;
        coef[4]=  0;
        coef[5]=  -1.0/2*CFL*e_;

        coef[6]= CFL*visc_*d_/dx;
        coef[7]= -2*CFL*visc_*d_/dx;
        coef[8]= CFL*visc_*d_/dx;

        break;
    case 3: // Finite Difference/Volume
        /*[1/8, 3/4, 1/8, 1/2*C, 0, -1/2*C]*/
        ic=0;
        coef[ic++]= 0;
        coef[ic++]= 0;
        coef[ic++]= 1;
        coef[ic++]= 0;
        coef[ic++]= 0;

        coef[ic++]=  - CFL*e_/12;
        coef[ic++]=  2*CFL*e_/3;
        coef[ic++]=  0;
        coef[ic++]=  -2*CFL*e_/3;
        coef[ic++]=  CFL*e_/12;

        coef[ic++]= - CFL*visc_*d_/dx/12;
        coef[ic++]= 4*CFL*visc_*d_/dx/3;
        coef[ic++]= -5*CFL*visc_*d_/dx/2;
        coef[ic++]= 4*CFL*visc_*d_/dx/3;
        coef[ic++]=  - CFL*visc_*d_/dx/12;
	break;
    case 4: //LS-FEM -  Peicewise parabolic
        // [1/10, 4/5, 1/10, 0, 0,]
        coef[0] = 0.1;
        coef[1] = 0.8;
        coef[2] = 0.1;
        coef[3] = 0.0;
        coef[4] = 0.0;

//    1/2*C, 0, -1/2*C, 0, 0
//	    ?negate
        coef[5] = 0.5*CFL;
        coef[6] = 0.0;
        coef[7] = -0.5*CFL;
        coef[8] = 0.0;
        coef[9] = 0.0;

//[   [-1/10, 1/5, 4/5, 1/5, -1/10,

        coef[10] = -0.1;
        coef[11] = 0.2;
        coef[12] = 0.8;
        coef[13] = 0.2;
        coef[14] = -0.1;
//
//     -1/4*C, C, 0, -C, 1/4*C]
//	    ?negate
        coef[15] = -0.25*CFL;
        coef[16] = CFL;
        coef[17] = 0;
        coef[18] = -CFL;
        coef[19] = 0.25*CFL;

        break;
    case 5: //LS-FEM -  Fourier 4
        double p2,p3,den0,den1,den2,den3;
        int ic;
        p2=pi*pi;
        p3=p2*pi;
        den0=(4.0/p2 + 1.0);
        den1=(30.0/p2 + 1.0);
        den2=(14.0/p2 + 3.0);
        den3=(32.0/p2 + 9.0);
        ic=0;
        coef[ic++] = 1.0/6.0/den0;
        coef[ic++] = 0.0;
        coef[ic++] = 4/(den0*p3);
        coef[ic++] = 0.0;
        coef[ic++] = 2.0/3.0/den0;
        coef[ic++] = 4/(den0*p2);
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = 1.0/6.0/den0;
        coef[ic++] = 0.0;
        coef[ic++] = -4/(den0*p3);
        coef[ic++] = 0.0;

        coef[ic++] = 0.25*CFL/den0;
        coef[ic++] = CFL/(den0*p2);
        coef[ic++] = 0.5*CFL/(den0*pi);
        coef[ic++] = 0.25*CFL/(den0*pi);
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = -CFL/(den0*pi);
        coef[ic++] = -0.5*CFL/(den0*pi);
        coef[ic++] = -0.25*CFL/den0;
        coef[ic++] = -CFL/(den0*p2);
        coef[ic++] = 0.5*CFL/(den0*pi);
        coef[ic++] = 0.25*CFL/(den0*pi);


        coef[ic++] = 0.0;
        coef[ic++] = 0.5*(3/p2 - 1)/den1;
        coef[ic++] = 0.0;
        coef[ic++] = -112.0/9.0/(den1*p3);
        coef[ic++] = 24/(den1*p2);
        coef[ic++] = (3/p2 + 2)/den1;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = 0.5*(3/p2 - 1)/den1;
        coef[ic++] = 0.0;
        coef[ic++] = 112.0/9.0/(den1*p3);

        coef[ic++] = 6*CFL/(den1*p2);
        coef[ic++] = -0.75*CFL/den1;
        coef[ic++] = 0.25*pi*CFL/den1;
        coef[ic++] = -2*CFL/(den1*pi);
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = -pi*CFL/den1;
        coef[ic++] = -28.0/3.0*CFL/(den1*pi);
        coef[ic++] = -6*CFL/(den1*p2);
        coef[ic++] = 0.75*CFL/den1;
        coef[ic++] = 0.25*pi*CFL/den1;
        coef[ic++] = -2*CFL/(den1*pi);


        coef[ic++] = -72/(den2*p3);
        coef[ic++] = 0.0;
        coef[ic++] = -1.5*(3/p2 + 1)/den2;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = -3*(3/p2 - 2)/den2;
        coef[ic++] = 32/(den2*p2);
        coef[ic++] = 72/(den2*p3);
        coef[ic++] = 0.0;
        coef[ic++] = -1.5*(3/p2 + 1)/den2;
        coef[ic++] = 0.0;

        coef[ic++] = -9*CFL/(den2*pi);
        coef[ic++] = -0.75*pi*CFL/den2;
        coef[ic++] = -2.25*CFL/den2;
        coef[ic++] = -80.0/3.0*CFL/(den2*p2);
        coef[ic++] = 18*CFL/(den2*pi);
        coef[ic++] = 3*pi*CFL/den2;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = -9*CFL/(den2*pi);
        coef[ic++] = -0.75*pi*CFL/den2;
        coef[ic++] = 2.25*CFL/den2;
        coef[ic++] = 80.0/3.0*CFL/(den2*p2);


        coef[ic++] = 0.0;
        coef[ic++] = 112.0/3.0/(den3*p3);
        coef[ic++] = 0.0;
        coef[ic++] = 0.375*(3/p2 + 4)/den3;
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = 32/(den3*p2);
        coef[ic++] = -0.75*(3/p2 - 8)/den3;
        coef[ic++] = 0.0;
        coef[ic++] = -112.0/3.0/(den3*p3);
        coef[ic++] = 0.0;
        coef[ic++] = 0.375*(3/p2 + 4)/den3;

        coef[ic++] = -4.5*CFL/(den3*pi);
        coef[ic++] = 6*CFL/(den3*pi);
        coef[ic++] = -80.0/3.0*CFL/(den3*p2);
        coef[ic++] = 2.25*CFL/den3;
        coef[ic++] = 9*CFL/(den3*pi);
        coef[ic++] = 28*CFL/(den3*pi);
        coef[ic++] = 0.0;
        coef[ic++] = 0.0;
        coef[ic++] = -4.5*CFL/(den3*pi);
        coef[ic++] = 6*CFL/(den3*pi);
        coef[ic++] = 80.0/3.0*CFL/(den3*p2);
        coef[ic++] = -2.25*CFL/den3;


    }
    std::cout << "CFL = " << CFL << std::endl;
    std::cout << " Coefficients\n";
    for(int i= 0; i<ncoef; i++) {
        std::cout << "  "<< coef[i] ;
    }
    std::cout << std::endl;

}

void RKWidget::setCFL(const double value) {
    if(CFL == value) return;
    dirty = true;
    CFL = value;
    //updateCoef(method);
}

void RKWidget::fillA() {
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

void RKWidget::blockFillA() {
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
                temp = coef[ccol+nb*nvar*2];
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
                    temp = coef[ccol+nb*nvar*2];
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

void RKWidget::fillB(int stg) {
    int ui;

    for( int i =0; i<N_; i += nblock ) {
        for(int nb = 0; nb < nblock; nb++) rhsb[i+nb] = U_[i];
    }
    for(int s=0 ; s<stg; s++) {
        if( b_a[nStage*stg+s] != 0.0) {
            for( int i =0; i<N_; i += nblock ) {
                for(int nb = 0; nb < nblock; nb++) rhsb[i+nb] += b_a[nStage*stg+s]*b_k[s*N_+i+nb];
            }
        }
    }
    Efunc(rhsb);
    Dfunc(rhsb);


    for( int i =0; i<N_; i += nblock ) {
        for(int nb = 0; nb < nblock; nb++) {
            rhsb[i+nb] = 0;
            ui = i-ndn;
            if(ui < 0) ui = N_+ui;
            //std::cout << i << " U_[n]  \n";
            for( int j =0; j<nvar; j++ ) {
                //std::cout << ui << "  ";
                rhsb[i+nb] += coef[nb*nvar*2 + nvar + j]*e_*E_[ui];
                rhsb[i+nb] += coef[nb*nvar*2 + 2*nvar + j]*d_*D_[ui];
                //std::cout << "coef [ "<< nb*nvar*2 + nvar + j << " ] * " << ui << std::endl;
                ui++;
                if(ui == N_) ui=0;
            }
        }
    }
}


void RKWidget::setMethod(int index) {
    if( index == method ) return;
    method = index;
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
    switch(method) {
    case 0: // Mid-Point
//	y1 = y0 + h f( y0+h/2 f(y0))
//	nStage = 2
//	b_a = matrix(SR,[[0,0],[1/2,0]])
//	b_b=[0,1]
//	b_c=[0,1/2]
        setTitle( tr("Mid-Point"));
        setNStage(2);
        b_a[nStage]=0.5;
        b_b[1]=1.0;
        b_c[1]=0.5;
        break;
    case 1: // Trapezoid
        //y1=y0+h/2f(y0) + h/2 f(y0+h f(y0))
//	nStage = 2
//	b_a = matrix(SR,[[0,0,],[1,0])
//	b_b = [1/2,1/2]
//	b_c=[0,1]
        setTitle( tr("Trapezoid"));
        setNStage(2);
        b_a[nStage]=1.0;
        b_b[0]=0.5;
        b_b[1]=0.5;
        b_c[1]=1.0;

        break;
    case 2: // RK4
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
        setTitle( tr("RK4"));
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
    case 3: // Arb-RK-4
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
        setTitle( tr("Arb-RK-4"));
        setNStage(4);
        b_a[4]=0.25;
        b_a[9]=1.0/3.0;
        b_a[14]=0.5;
        b_b[3]=1.0;
        b_c[1]=0.25;
        b_c[2]=1.0/3.0;
        b_c[3]=0.5;
        break;
    default:
        setMethod(0);
    }

}

void RKWidget::setNStage(int arg1) {
    if(nStage != arg1) {
        if(nStage != 0) {
            delete[] b_a;
            delete[] b_b;
            delete[] b_c;
            if( N_ != 0 ) delete[] b_k;
        }
        b_a = new double[arg1*arg1];
        b_b = new double[arg1];
        b_c = new double[arg1];
        if( N_ != 0) b_k = new double[arg1*N_];
        nStage = arg1;
    }
    for(int i=0; i<nStage; i++ ) {
        for(int j=0; j<nStage; j++) b_a[j*nStage+i] = 0.0;
        b_b[i] = 0.0;
        b_c[i] = 0.0;
    }
}
bool RKWidget::canSolve(int equ) {
    return true;
}





