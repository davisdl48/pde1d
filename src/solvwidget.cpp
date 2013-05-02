/*
    Base Class for solvers - provides the interface to pde1d
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


#include "solvwidget.h"
#include <iostream>
#include <cmath>
#include <qwt_symbol.h>

SolvWidget::SolvWidget ( QWidget *parent ) :  QDockWidget ( parent )
{
    setupUi ( );
    curve = new QwtPlotCurve();
    curve->setRenderHint(QwtPlotItem::RenderAntialiased, true);
    curve->setLegendAttribute(QwtPlotCurve::LegendShowBrush, true);
    curve->setLegendAttribute( QwtPlotCurve::LegendShowLine , true);
    curve->setLegendAttribute( QwtPlotCurve::LegendShowSymbol, true);
    connect ( plotColor, SIGNAL ( valueChanged ( QColor ) ), this, SLOT ( setColor ( QColor ) ) );
    connect ( plotNameEdit, SIGNAL ( textEdited ( QString ) ), this, SLOT ( setTitle ( QString ) ) );
    dataNames.append(tr("X")); // locations
    dataNames.append(tr("U")); // default dependant variable
    dataNames.append(tr("E")); // convective flux - e.g. c*U or 1/2 U^2
    dataNames.append(tr("D")); // diffusive flux - e.g. nu*U
    dataNames.append(tr("Jac")); // Jacobian d(D_x)/dU
    N_ = 0;
    dt = 0.0;
    pi = 4 * atan ( 1.0 );
    cStep = 0;
    dirty = true;
    setFeatures ( QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetClosable);
    samset = false;
    unstable = false;
    nghost = 0;
    nfghost = 1;
}

SolvWidget::SolvWidget ( const SolvWidget& other )
{

}

SolvWidget::~SolvWidget()
{
    QHash<QString, double *>::const_iterator iter = data.constBegin();
    while( iter != data.constEnd()) {
      delete[] iter.value();
      iter++;
    }
    data.clear();
    curve->detach();
}

SolvWidget& SolvWidget::operator= ( const SolvWidget& other )
{
    return *this;
}

bool SolvWidget::operator== ( const SolvWidget& other ) const
{
    return (this == &other);
}


QColor SolvWidget::getPenColor()
{
    return plotColor->getValue();
}

QString& SolvWidget::getTitle()
{
    if(unstable) {
      QString utitle = "*"+title+"*";
      return utitle;
    }
    return title;
}

void SolvWidget::setColor ( QColor color )
{
    const QwtSymbol *sym;
    QwtSymbol::Style style;
    QBrush brush;
    QPen pen;
    QSize size;
    plotColor->setValue ( color );
    pen = curve->pen();
    pen.setColor(color);
    curve->setPen(pen);
    sym = curve->symbol();
    if( sym != NULL && sym->style() != QwtSymbol::NoSymbol ) {
      style = sym->style();
      brush = sym->brush();
      brush.setColor(color);
      pen = sym->pen();
      pen.setColor(color);
      size = sym->size();
      curve->setSymbol( new QwtSymbol(style,brush,pen,size) );    
    }
}

void SolvWidget::setTitle ( QString newTitle )
{
    title = newTitle;
    plotNameEdit->setText ( title );
    setWindowTitle(title);
}

void SolvWidget::initSin ( const double value )
{
    samset = false;
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

void SolvWidget::setCFL ( const double value )
{
    CFL = value;
    if ( dx != 0.0 ) {
        dt = CFL * dx;
    }
    dirty=true;
}

void SolvWidget::setSpeed ( const double value )
{
    double c;
    c = value;
}

void SolvWidget::setup ( const size_t size, const double cycles, const double cfl )
{
    this->cycles = cycles;
    setSize ( size );
    //initSin ( cycles );
    setCFL ( cfl );
}

const size_t SolvWidget::getSize()
{
    return N_;
}

const double SolvWidget::getCFL()
{
    return CFL;
}


double* SolvWidget::getX()
{
    return X_;
}

double* SolvWidget::getU()
{
    return U_;
}

int SolvWidget::getCurrentStep()
{
    return cStep;
}

double* SolvWidget::getIdeal()
{
    double travel;
    dx = 2 * pi / N_;
    double amp;
    if(e_ == 0.0) {
      travel = 0.5;
    }else{
      travel = totCFL/ ( double ) N_;
    }
    travel = travel - std::floor ( travel );
    amp = exp(-(totCFL-N_/2.0)*dx*visc_*d_*cycles*cycles);
    double vx0 = ( 1 - travel ) * 2 * pi;
    double vx;
    for ( size_t i = 0; i < N_; i++ ) {
        vx = vx0 + X_[i];
        if ( vx >= 2 * pi ) vx -= ( 2 * pi );
        Ideal_[i] = amp*sin ( cycles * vx );
    }
    return Ideal_;
}

double SolvWidget::getCycles()
{
    return cycles;
}

double SolvWidget::getTravel()
{
    return totCFL / ( double ) N_;
}

void SolvWidget::closeEvent ( QCloseEvent* ev )
{
    emit dockClose(id);
    ev->accept();
}

int SolvWidget::getId()
{
    return id;
}

void SolvWidget::setId ( const int value )
{
    id = value;
}

void SolvWidget::Efunc(double* Udat) {
   switch( equation ) {
     default:
     case 0:
     case 1:
        for(size_t nn=0; nn<N_; nn++) {
            E_[nn] = Udat[nn];
	    J_[nn] = 1.0;
        }
        return;
     case 2:
     case 3:
        for(size_t nn=0; nn<N_; nn++) {
            E_[nn] = Udat[nn]*Udat[nn]/2.0;
	    J_[nn] = Udat[nn];
        }
        return;
   }
}

void SolvWidget::Dfunc(double* Ddat) {
    for(size_t nn=0; nn<N_; nn++) {
        D_[nn]= Ddat[nn];
    }
}


void SolvWidget::setEquation(int index) {
    equation = index;
    dirty=true;
    switch(index) {
    default:
    case 0:
        d_=0.0;
        e_=1.0;
        break;
    case 1:
        d_=1.0;
        e_=0.0;
        break;
    case 2:
        e_=1.0;
        d_=0.0;
        break;
    case 3:
        e_=1.0;
        d_=1.0;
    }
}

void SolvWidget::setViscosity(double value) {
  if(value == visc_) return;
    visc_=value;
    dirty=true;   
}

void SolvWidget::resize(int value) {
    if ( value == N_ ) return;
    cStep = 0;
    QHash<QString, double *>::const_iterator iter = data.constBegin();
    while( iter != data.constEnd()) {
      delete[] iter.value();    
      iter++;
    }
    QStringList::const_iterator constIterator;
    for (constIterator = dataNames.constBegin(); constIterator != dataNames.constEnd();
            ++constIterator){
         data.insert(*constIterator, new double[value]);
    }
    
    QHash<QString, double *>::const_iterator fiter = fluxes.constBegin();
    while( fiter != fluxes.constEnd()) {
      delete[] fiter.value();    
      fiter++;
    }
    QStringList::const_iterator constfIterator;
    for (constfIterator = fluxNames.constBegin(); constfIterator != fluxNames.constEnd();
            ++constfIterator){
         fluxes.insert(*constfIterator, new double[value+1]);
    }
    
  // temporary for compatablitiy  
    if ( N_ != 0 ) {
	delete[] Ideal_;
    }
    N_ = value;
    Ideal_ = new double[N_];
    
    U_ = data.value(tr("U"),0);
    X_ = data.value(tr("X"),0);
    E_ = data.value(tr("E"),0);
    J_ = data.value(tr("Jac"),0);
    D_ = data.value(tr("D"),0);
    
    // End temporary data values 
    initSin ( cycles );
}

void SolvWidget::setupUi() {
    if (objectName().isEmpty())
        setObjectName(QString::fromUtf8("SolvWidget"));
    dockWidgetContents = new QWidget();
    dockWidgetContents->setObjectName(QString::fromUtf8("dockWidgetContents"));
    verticalLayout = new QVBoxLayout(dockWidgetContents);
    verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
    plotNameLabel = new QLabel(dockWidgetContents);
    plotNameLabel->setText(tr("Title"));
    plotNameLabel->setObjectName(QString::fromUtf8("plotNameLabel"));

    verticalLayout->addWidget(plotNameLabel);

    plotNameEdit = new QLineEdit(dockWidgetContents);
    plotNameEdit->setMinimumWidth(50);
    plotNameEdit->setMaximumWidth(222);
    plotNameEdit->setObjectName(QString::fromUtf8("plotNameEdit"));

    verticalLayout->addWidget(plotNameEdit);

    plotColorLabel = new QLabel(dockWidgetContents);
    plotColorLabel->setText(tr("Plot Color"));
    plotColorLabel->setObjectName(QString::fromUtf8("plotColorLabel"));

    verticalLayout->addWidget(plotColorLabel);

    plotColor = new MyColorButton(dockWidgetContents);
    plotColor->setMinimumWidth(30);
    plotColor->setMaximumWidth(100);
    plotColor->setObjectName(QString::fromUtf8("plotColor"));

    verticalLayout->addWidget(plotColor);

    verticalSpacer = new QSpacerItem(240, 175, QSizePolicy::Minimum, QSizePolicy::Expanding);

    verticalLayout->addItem(verticalSpacer); 

    setWidget(dockWidgetContents);
    
    //dockWidgetContents->resize(250, 640);
    
    //dockWidgetContents->setMaximumWidth ( 250 );

    //retranslateUi(SolvWidget);

    //QMetaObject::connectSlotsByName(SolvWidget);
}

bool SolvWidget::canSolve(int equ) {
    return false;
}

bool SolvWidget::isUnstable() {
    return unstable;
}

bool SolvWidget::isOK() {
    return (!unstable)&&canSolve(equation);
}

QwtPlotCurve* SolvWidget::getCurve(QString xvalue, QString value ) {
    double xx[1] = {0.0};
    double uu[1] = {0.0};
    if(unstable) {
      curve->setTitle( "*"+title+"*" );
      if( !samset ) curve->setSamples ( xx,uu,1);
      //curve->setSamples ( xx,uu,1);
    }else{
      curve->setTitle( title );
      if(!id) getU();
      if( !samset ) curve->setSamples ( data.value(xvalue,X_), data.value(value,U_), N_ );
      //curve->setSamples ( data.value(xvalue,X_), data.value(value,U_), N_ );
    }
    samset = true;
    return curve;
}

void SolvWidget::setSize(const size_t size) {
    resize(size);
}
