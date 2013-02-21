/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2012  Davis Family <email>

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

SolvWidget::SolvWidget ( QWidget *parent ) :  QDockWidget ( parent ), Ui_SolvWidget()
{
    setupUi ( this );
    connect ( plotColor, SIGNAL ( activated ( QColor ) ), this, SLOT ( setColor ( QColor ) ) ) ;
    connect ( plotNameEdit, SIGNAL ( textEdited ( QString ) ), this, SLOT ( setTitle ( QString ) ) );
    N = 0;
    dt = 0.0;
    pi = 4 * atan ( 1.0 );
    cStep = 0;
    frameNum = 0;
    setFeatures ( QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetClosable);
}

SolvWidget::SolvWidget ( const SolvWidget& other )
{

}

SolvWidget::~SolvWidget()
{
    std::cout << "Delete solve widget\n";
    if ( N != 0 ) {
        delete[] U;
        delete[] ideal;
        delete[] x;
    }
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
    return col;
}

QString& SolvWidget::getTitle()
{
    return title;
}

void SolvWidget::setColor ( QColor color )
{
    col = color;
    plotColor->setColor ( col );
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
    if ( N == 0 ) setSize ( 100 );
    totCFL = N / 2.0;
    cycles = value;
    dx = 2 * pi / N;
    for ( size_t i = 0; i < N; i++ ) {
        x[i] = dx * i;
        //U[i] = sin ( cycles * x[i] + travel);
    }
    getIdeal();
    for ( size_t i = 0; i < N; i++ ) {
        U[i] = ideal[i];
    }

}

void SolvWidget::setCFL ( const double value )
{
    CFL = value;
    if ( dx != 0.0 ) {
        dt = CFL * dx;
    }
}

void SolvWidget::setSpeed ( const double value )
{
    double c;
    c = value;
}

void SolvWidget::setup ( const size_t size, const double cycles, const double cfl )
{
    setSize ( size );
    initSin ( cycles );
    setCFL ( cfl );
}

const size_t SolvWidget::getSize()
{
    return N;
}

const double SolvWidget::getCFL()
{
    return CFL;
}


double* SolvWidget::getX()
{
    return x;
}

double* SolvWidget::getU()
{
    return U;
}

int SolvWidget::getCurrentStep()
{
    return cStep;
}

double* SolvWidget::getIdeal()
{
    double travel;
    dx = 2 * pi / N;
    travel = totCFL / ( double ) N;
    travel = travel - std::floor ( travel );
    double vx0 = ( 1 - travel ) * 2 * pi;
    double vx;
    for ( size_t i = 0; i < N; i++ ) {
        vx = vx0 + x[i];
        if ( vx >= 2 * pi ) vx -= ( 2 * pi );
        ideal[i] = sin ( cycles * vx );
    }
    return ideal;
}

double SolvWidget::getCycles()
{
    return cycles;
}

double SolvWidget::getTravel()
{
    return totCFL / ( double ) N;
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



