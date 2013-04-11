/*
    A QDockWidget to hold a table for changing plot curve attributes
    such as line style, symbol, color.
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


#ifndef CURVETABDOCK_H
#define CURVETABDOCK_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDockWidget>
#include <QtGui/QHeaderView>
#include <QtGui/QTableWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include <QLabel>
#include <QList>
#include "solvwidget.h"

class CurvesModel;

class CurveTabDock : public QDockWidget
{
    Q_OBJECT

public:
    explicit CurveTabDock(const QString& title, QWidget* parent = 0, Qt::WindowFlags flags = 0);
    explicit CurveTabDock(QWidget* parent = 0, Qt::WindowFlags flags = 0);
    virtual ~CurveTabDock();
    
    void setCurvesModel(CurvesModel * newModel) ;




protected:
  void setupUi(const CurveTabDock* ctd = NULL) ;
        QWidget *dockWidgetContents;
        QHBoxLayout *verticalLayout;
        QTableView *curveTable;
	CurvesModel *curveModel;
	
	
    };

#endif // CURVETABDOCK_H
