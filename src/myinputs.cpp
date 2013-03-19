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


#include "myinputs.h"
#include <qcolordialog.h>

MyIntInput::MyIntInput(QWidget* parent): QLineEdit(parent )
{
    val = new QIntValidator(1,9999999,this);
    this->setValidator(val);
    if( text().isEmpty() ) {
        myValue=1;
    } else {
        myValue=text().toInt();
    }
    setText( QString("%1").arg(myValue));
    connect(this,SIGNAL(editingFinished()),this,SLOT(done()));
}

MyIntInput::MyIntInput(const QString& value, QWidget* parent): QLineEdit(value,parent )
{
    val = new QIntValidator(1,9999999,this);
    this->setValidator(val);
    if( text().isEmpty() ) {
        myValue=1;
    } else {
        myValue=text().toInt();
    }
    setText( QString("%1").arg(myValue));
    connect(this,SIGNAL(editingFinished()),this,SLOT(done()));

}

int MyIntInput::getValue() {
    return myValue;
}

void MyIntInput::done() {
    bool ok;
    myValue = text().toInt(&ok);
    if(ok) {
        emit valueChanged(myValue);
    } else {
        setText( QString("%1").arg(myValue));
    }
}

void MyIntInput::setValue(int value) {
    myValue = value;
    setText( QString("%1").arg(myValue));
}


MyDoubInput::MyDoubInput(QWidget* parent): QLineEdit(parent )
{
    val = new QDoubleValidator(this);
    this->setValidator(val);
    if( text().isEmpty() ) {
        myValue=0.0;
    } else {
        myValue=text().toDouble();
    }
    setText( QString("%1").arg(myValue));
    connect(this,SIGNAL(editingFinished()),this,SLOT(done()));
}

MyDoubInput::MyDoubInput(const QString& value, QWidget* parent): QLineEdit(value,parent)
{
    val = new QDoubleValidator(this);
    this->setValidator(val);
    if( text().isEmpty() ) {
        myValue=0.0;
    } else {
        myValue=text().toDouble();
    }
    setText( QString("%1").arg(myValue));
    connect(this,SIGNAL(editingFinished()),this,SLOT(done()));
}



void MyDoubInput::done() {
    bool ok;
    myValue = text().toDouble(&ok);
    if(ok) {
        emit valueChanged(myValue);
    } else {
        setText( QString("%1").arg(myValue));
    }
}

void MyDoubInput::setValue(double value) {
    myValue = value;
    setText(QString("%1").arg(myValue));
}

MyColorButton::MyColorButton(QWidget* parent) : QPushButton(parent)
{
        myValue = Qt::red;
        setPalette(QPalette(myValue));
        setAutoFillBackground(true);
	connect(this,SIGNAL(clicked()),this,SLOT(setColor()));
}

MyColorButton::MyColorButton(const QColor& value, QWidget* parent) : QPushButton(parent)
{
  if (value.isValid()) {
        //setText(color.name());
        myValue=value;
    }else{
      myValue = Qt::red;
    }
        setPalette(QPalette(myValue));
        setAutoFillBackground(true);
	connect(this,SIGNAL(clicked()),this,SLOT(setColor()));
}

QColor MyColorButton::getValue() {
    return myValue;
}

void MyColorButton::setValue(QColor value) {
    if (value.isValid()) {
        //setText(color.name());
        setPalette(QPalette(value));
        setAutoFillBackground(true);
        myValue=value;
    }
}

void MyColorButton::setColor() {
    bool native;
    native = true;
    QColor color;
    QWidget *qw;
    qw = this;
    if (native) {
        color = QColorDialog::getColor(myValue, qw, tr("Select Color"), QColorDialog::ShowAlphaChannel);
    } else {
        color = QColorDialog::getColor(myValue, qw, tr("Select Color"), QColorDialog::DontUseNativeDialog | QColorDialog::ShowAlphaChannel);
    }

    if (color.isValid()) {
        //setText(color.name());
        setPalette(QPalette(color));
        setAutoFillBackground(true);
	emit valueChanged(color);
    }
}
