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


#ifndef ERRTABDOCK_H
#define ERRTABDOCK_H

#include <QtGui/qdockwidget.h>
#include <ui_errtabdock.h>


class ErrTabDock : public QDockWidget, public Ui_ErrTabDock
{
  Q_OBJECT

public:
    explicit ErrTabDock ( const QString& title, QWidget* parent = 0, Qt::WindowFlags flags = Qt::Widget );
    explicit ErrTabDock ( QWidget* parent = 0, Qt::WindowFlags flags = Qt::Widget );
    void errTab();
};

#endif // ERRTABDOCK_H
