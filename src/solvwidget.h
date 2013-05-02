/*
    This is a base class for all solution methods.
    For new methods override the constructor to provide any editing widgets,
    setSize(size) to allocate custom data, the destructor to delete custom data
    and step(nsteps) to solve.
    Copyright (C) 2012  Davis Family davisdl48@gmail.com

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


#ifndef SOLVWIDGET_H
#define SOLVWIDGET_H
#include <QtGui/qcolor.h>
#include <QtGui/QCloseEvent>
#include <QtGui/QDockWidget>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include <qwt_plot_curve.h>
#include "myinputs.h"

/// Base widget to be inherited to add a numerical solver
/**
This base class provides much of the structure for interaction
with the remainder of the code.  To add a new solver inherit from this class
and override the functions.

\code
void step(const size_t nStep ) {  -- to solve for nStep times }
bool canSolve(int equ); // indicates the solver can solve equation number equ
\endcode

The constructor should add variable names to dataNames for data at \f$ x_i \f$ nodes
and fluxNames for data at midpoints, \f$ x_i-1/2 \f$.  e.g.
\code
dataNames.append(tr("MyVariableName"));
\endcode
This will make the variable double pointer available as
\code
double *myVariable;
myVariable = data.value(tr("MyVariableName"));
\endcode

\todo make all variables accessible for plotting by name

default names include 
 "X" // locations
 "U" // default dependant variable
 "E" // convective flux - e.g. c*U or 1/2 U^2
 "D" // diffusive flux - e.g. nu*U
 "Jac" // Jacobian dE/dU
  
To add custom widgets/controls for the solver, create QWidget items, connect the signals
to slots and add them to verticalLayout within the constructor.
\code
    upwLabel = new QLabel(tr("Upwinding"));
    upwInput = new MyDoubInput(0.0,this,-0.1,1.1,0.01,6);// (Initial Value, this, minValue, maxValue, increment, precision)
    connect(upwInput,SIGNAL(valueChanged(double)),this,SLOT(setUpwind(double)));
    verticalLayout->addWidget(upwLabel);
    verticalLayout->addWidget(upwInput);
    .
    .
    .
    verticalLayout->addItem(verticalSpacer);
\endcode

Also, initialize any custom parameters in the constructor and
set the initial values in any custom widgets.

For additional data that does not fit the form as a double array (double pointer)
override setSize(int value) to allocate the custom data items. 
\note at the end of the setSize() method make sure to call resize(int value) to
perform the automatic size updates.  resize() will allocate all named variables in dataNames
and fluxNames and reset the array size variable N_.

Override the destructor to delete any custom data types.

The new solver must currently be added to pde1d.cpp
Include the new headers
\code
#include "specwidget.h"
#include "pswidget.h"
\endcode
in the constructor add the new solvers to the combobox, e.g.
\code
    solvModel->appendRow( new QStandardItem ( tr("Spectral FFT") ));    // addSolver(8)
    solvModel->appendRow( new QStandardItem ( tr("Pseudo Spectral") )); // addSolver(9)
\endcode
and in addSolver ( int index ) add the additionnal case, e.g.
\code    
    case 8:
        SpecWidget *specw;
        specw= new SpecWidget ( this );
        addIt(specw);
        break;
    case 9:
        PSWidget *psw;
        psw = new PSWidget ( this );
        addIt(psw);
        break;
\endcode

\todo make the new solver a shared object and load it with dlopen so pde1d does not need to be modified.
*/
class SolvWidget : public  QDockWidget
{
    Q_OBJECT

public:
    SolvWidget ( QWidget *parent = 0 );
    SolvWidget ( const SolvWidget& other );
    virtual ~SolvWidget();
    virtual SolvWidget& operator= ( const SolvWidget& other );
    virtual bool operator== ( const SolvWidget& other ) const;
    QString& getTitle();
    
    QwtPlotCurve * getCurve( QString xvalue= tr("X"), QString value = tr("U") );
      
    QColor getPenColor();

    const size_t getSize() ;
    void resize(int value );
    
    virtual void initSin ( const double cycles = 1.0 );

    virtual void setCFL ( const double value = 1.0 );
    const double getCFL();
    void setSpeed ( const double value = 1.0 );
    const double getSpeed();

    double * getX();
    virtual double * getU();
    double * getIdeal();

    void setup ( const size_t size = 100, const double cycles = 1.0, const double cfl = 1.0 );
    int getCurrentStep();
    double getCycles();
    double getTravel();
    

    virtual bool canSolve(int equ);
    bool isUnstable() ;
    bool isOK();

    int getId() ;
    void setId (const int value ) ;

    void closeEvent ( QCloseEvent *ev );
    
    virtual void setSize ( const size_t size = 100 );
    virtual void step ( const size_t nStep ) = 0;
    
    QWidget *dockWidgetContents;
    QVBoxLayout *verticalLayout;
    QLabel *plotNameLabel;
    QLineEdit *plotNameEdit;
    QLabel *plotColorLabel;
    MyColorButton *plotColor;
    QSpacerItem *verticalSpacer;
    void setupUi();
    QwtPlotCurve *curve;
   
    virtual void setEquation( int index) ;

public slots:
    void setColor ( QColor color ) ;
    void setTitle ( QString title ) ;
    
    void setViscosity( double value = 0.0 ) ;
    

signals:
  void dockClose(int id);


protected:
    QString title;
    QHash<QString,double *> data;
    QStringList dataNames;
    QHash<QString,double *> fluxes;
    QStringList fluxNames;

    size_t N_;
    double *U_;
    double *X_;
    /// U_t + e_*E_(U)_x  + d_*D_(U)_xx = 0.0
    double *E_; // Flux Vector- scalar
    double *J_; // Jacobian dE/dU
    double *D_; // Dissipation Variable
    double *Ideal_;
    double d_; // dissipation coefficient = viscosity*dx/c*CFL
    double e_; // flux coefficient = CFL
    void Efunc(double *Udat);
    void Dfunc(double *Ddat);
    double visc_; // viscosity
    double speed_;
    double dt;
    double dx;
    double CFL;
    double cycles;
    int cStep;
    double pi;
    double totCFL;
    int id;
    bool dirty;
    int equation;
    bool unstable;
    bool samset;
    int nghost;
    int nfghost;
};

#endif // EULEREXWIDGET_H
