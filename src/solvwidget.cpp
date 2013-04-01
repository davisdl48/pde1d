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

SolvWidget::SolvWidget ( QWidget *parent ) :  QDockWidget ( parent )
{
    setupUi ( );
    connect ( plotColor, SIGNAL ( valueChanged ( QColor ) ), this, SLOT ( setColor ( QColor ) ) ) ;
    connect ( plotNameEdit, SIGNAL ( textEdited ( QString ) ), this, SLOT ( setTitle ( QString ) ) );
    dataNames.append(tr("X"));
    dataNames.append(tr("U"));
    dataNames.append(tr("E")); // convective flux - e.g. c*U or 1/2 U^2
    dataNames.append(tr("D")); // diffusive flux - e.g. nu*U
    dataNames.append(tr("Jac")); // Jacobian d(D_x)/dU
    N_ = 0;
    dt = 0.0;
    pi = 4 * atan ( 1.0 );
    cStep = 0;
    frameNum = 0;
    dirty = true;
    setFeatures ( QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetClosable);
    
    unstable = false;
}

SolvWidget::SolvWidget ( const SolvWidget& other )
{

}

SolvWidget::~SolvWidget()
{
    //std::cout << "Delete solve widget\n";
    QHash<QString, double *>::const_iterator iter = data.constBegin();
    while( iter != data.constEnd()) {
      std::cout << iter.key().toLocal8Bit().constData() << ":  " << iter.value() << std::endl;
      delete[] iter.value();
      iter++;
    }
    data.clear();
}

SolvWidget& SolvWidget::operator= ( const SolvWidget& other )
{
    return *this;
}

bool SolvWidget::operator== ( const SolvWidget& other ) const
{
    return (this == &other);
}


QColor SolvWidget::getColor()
{
    return plotColor->getValue();
}

QString& SolvWidget::getTitle()
{
    return title;
}

void SolvWidget::setColor ( QColor color )
{
    plotColor->setValue ( color );
}

void SolvWidget::setTitle ( QString newTitle )
{
    std::cout << "SolvWidget::setTitle( " << title.toLocal8Bit().data() << " )\n";
    title = newTitle;
    plotNameEdit->setText ( title );
    setWindowTitle(title);
}

void SolvWidget::initSin ( const double value )
{
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
    //std::cout << "Recieved close event\n" << std::flush;
    emit dockClose(id);
    ev->accept();
}

int SolvWidget::getId()
{
    return id;
}

void SolvWidget::setId ( int value )
{
    id = value;
    std::cout << " dock " << title.toLocal8Bit().data() << "  id = " << id << std::endl;
}

void SolvWidget::Efunc(double* Udat) {
    if(burg) {
        for(size_t nn=0; nn<N_; nn++) {
            E_[nn] = Udat[nn]*Udat[nn]/2.0;
	    J_[nn] = Udat[nn];
        }
    } else {
        for(size_t nn=0; nn<N_; nn++) {
            E_[nn] = Udat[nn];
	    J_[nn] = 1.0;
        }
    }
}

void SolvWidget::Dfunc(double* Ddat) {
    for(size_t nn=0; nn<N_; nn++) {
        D_[nn]= Ddat[nn];
    }
}


void SolvWidget::setEquation(int index) {
    equation = index;
    burg=false;
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
        burg=true;
        e_=1.0;
        d_=0.0;
        break;
    case 3:
        burg=true;
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
    initSin ( cycles );
}

bool SolvWidget::getBurg() {
    return burg;
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
    plotNameEdit->setObjectName(QString::fromUtf8("plotNameEdit"));

    verticalLayout->addWidget(plotNameEdit);

    plotColorLabel = new QLabel(dockWidgetContents);
    plotColorLabel->setText(tr("Plot Color"));
    plotColorLabel->setObjectName(QString::fromUtf8("plotColorLabel"));

    verticalLayout->addWidget(plotColorLabel);

    plotColor = new MyColorButton(dockWidgetContents);
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
double SolvWidget::getLineWidth() {
    return lineWidth;
}
void SolvWidget::setLineWidth(double lw) {
    lineWidth=lw;
}
QwtPlotCurve* SolvWidget::getCurve(QString value, QString xvalue) {
    curve->setTitle( title );
    curve->setSamples ( data.value(value,U_), data.value(xvalue,X_), N_ );
    // tobe replaced - change pen through gui
    curve->setPen ( QPen ( QBrush ( plotColor->getValue() ), 1.5, Qt::DashLine ) );
}
