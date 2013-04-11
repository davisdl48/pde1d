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


#ifndef MYINPUTS_H
#define MYINPUTS_H

#include <QtGui/QLineEdit>
#include <qvalidator.h>
#include <qcombobox.h>
#include <qpushbutton.h>


class MyIntInput : public QLineEdit
{
Q_OBJECT
public:
    MyIntInput(QWidget* parent = 0);
    MyIntInput(const QString& value, QWidget* parent = 0);
    
    QIntValidator * val;
    int getValue();
    int myValue;
  public   slots:
    void done();
    void setValue(int value) ;
    
signals:
   void valueChanged(int value);
};


class MyDoubInput : public QLineEdit
{
Q_OBJECT
public:
    MyDoubInput(QWidget* parent = 0);
    MyDoubInput(const QString& , QWidget* parent = 0);
    MyDoubInput(double value, QWidget *parent=0,double lower=-1.7e+308, double upper=1.7e+308,  double singleStep=0.01, int precision=6);
      
    QDoubleValidator * val;
    double getValue();
  public  slots:
    void done();
    void setValue(double value) ;
    
signals:
   void valueChanged(double value);
private:
    double myValue;
    double minValue;
    double maxValue;
    double increment;
    int prec;
};

class MyColorButton : public QPushButton
{
Q_OBJECT
public:
    MyColorButton(QWidget* parent = 0);
    MyColorButton(const QColor& , QWidget* parent = 0);
     
    QColor getValue();
public  slots:
    void setColor();
    void setValue(QColor value);
    
signals:
   void valueChanged(QColor value);
private:
    QColor myValue;
};

#endif // MYINTINPUT_H
