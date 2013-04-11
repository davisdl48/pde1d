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

//#include "eulerexwidget.moc"
/*! Base widget should be inherited to add a numerical solver
 * 
 * To add a new solver inherit from this class and override
 * ~ - the destructor to delete data
 * setSize(int ) - dimension data and delete if nessecary 
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
    double getLineWidth();
    void setLineWidth(double lw);

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
    
    bool getBurg();
   

    int getId() ;
    void setId ( int value ) ;

    virtual void closeEvent ( QCloseEvent *ev ) ;
    virtual void setSize ( const size_t size = 100 ) = 0;
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
    double dt;
    double dx;
    double CFL;
    double cycles;
    int cStep;
    double pi;
    double totCFL;
    int frameNum;
    int id;
    bool burg;
    bool dirty;
    int equation;
    bool unstable;
    bool samset;
    double lineWidth;
};

#endif // EULEREXWIDGET_H
