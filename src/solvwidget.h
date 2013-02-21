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


#ifndef SOLVWIDGET_H
#define SOLVWIDGET_H
#include <ui_solvwidget.h>
#include <QtGui/qcolor.h>
#include <QtGui/QCloseEvent>
#include "kcolorcombo.h"
//#include "eulerexwidget.moc"

class SolvWidget : public  QDockWidget, public Ui_SolvWidget
{
    Q_OBJECT

public:
    SolvWidget ( QWidget *parent = 0 );
    SolvWidget ( const SolvWidget& other );
    virtual ~SolvWidget();
    virtual SolvWidget& operator= ( const SolvWidget& other );
    virtual bool operator== ( const SolvWidget& other ) const;
    QColor getColor();
    QString& getTitle();

    const size_t getSize() ;
    void initSin ( const double cycles = 1.0 );

    virtual void setCFL ( const double value = 1.0 );
    const double getCFL();
    void setSpeed ( const double value = 1.0 );
    const double getSpeed() ;

    double * getX();
    double * getU();
    double * getIdeal();

    void setup ( const size_t size = 100, const double cycles = 1.0, const double cfl = 1.0 );
    int getCurrentStep();
    double getCycles();
    double getTravel();

    int getId() ;
    void setId ( int value ) ;

    virtual void closeEvent ( QCloseEvent *ev ) ;
    virtual void setSize ( const size_t size = 100 ) = 0;
    virtual void step ( const size_t nStep ) = 0;

public slots:
    void setColor ( QColor color ) ;
    void setTitle ( QString title ) ;

signals:
  void dockClose(int id);


protected:
    QColor col;
    QString title;

    size_t N;
    double *U;
    double *ideal;
    double *x;
    double dt;
    double dx;
    double CFL;
    double cycles;
    int cStep;
    double pi;
    double totCFL;
    int frameNum;
    int id;
};

#endif // EULEREXWIDGET_H
