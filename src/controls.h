/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2013  Davis Family <email>

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


#ifndef CONTROLS_H
#define CONTROLS_H

#include <QtGui/qdockwidget.h>
#include "ui_controls.h"


class Controls : public QDockWidget, public Ui_Controls
{
    Q_OBJECT

public:
    explicit Controls ( const QString& title = "Controls", QWidget* parent = 0, Qt::WindowFlags flags = 0 );
    //Controls();
    Controls ( const Controls& other );
    virtual ~Controls();
    virtual Controls& operator= ( const Controls& other );
    virtual bool operator== ( const Controls& other ) const;


};

#endif // CONTROLS_H
