/*
    DockWidget to hold a table to control curve plot attibutes.
    Copyright (C) 2013  Doug Davis  davisdl48@gmail.com

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


#include "curvetabdock.h"
#include "curvesmodel.h"
#include "itemdelegate.h"
#include "colordelegate.h"


CurveTabDock::CurveTabDock(const QString& title, QWidget* parent, Qt::WindowFlags flags) : QDockWidget(title, parent, flags)
{
    setObjectName(title);
  setupUi(this);
}

CurveTabDock::CurveTabDock(QWidget* parent, Qt::WindowFlags flags) : QDockWidget(parent, flags)
{
    setObjectName(QString::fromUtf8("Curve Properties"));
  setupUi(this);
}

CurveTabDock::~CurveTabDock() {
}

void CurveTabDock::setCurvesModel(CurvesModel *newModel) {
    curveTable->setModel(newModel);
    return;
}

void CurveTabDock::setupUi(const CurveTabDock* ctd) {
    dockWidgetContents = new QWidget();
    verticalLayout =  new QHBoxLayout(dockWidgetContents);
    curveTable = new QTableView(dockWidgetContents);
    //dockWidgetContents->setSizePolicy(QSizePolicy::Maximum,QSizePolicy::Maximum);
    
    curveTable->horizontalHeader()->setVisible(true);
    curveTable->verticalHeader()->setVisible(false);
    curveTable->setEditTriggers(QAbstractItemView::AllEditTriggers);
    
    ColorDelegate *penColor = new ColorDelegate(this);
    curveTable->setItemDelegateForColumn(1,penColor);
    curveTable->setColumnWidth(1,170);
    
    PenStyleDelegate *penStyle = new PenStyleDelegate(this);
    curveTable->setItemDelegateForColumn(2,penStyle);
    
    LineWidthDelegate *lineWidth = new LineWidthDelegate(this);
    curveTable->setItemDelegateForColumn(3,lineWidth);
    
    SymbolStyleDelegate *symStyle = new SymbolStyleDelegate(this);
    curveTable->setItemDelegateForColumn(4,symStyle);
    
    ColorDelegate *symColor = new ColorDelegate(this);
    curveTable->setItemDelegateForColumn(5,symColor);
    curveTable->setColumnWidth(5,170);
    
    verticalLayout->addWidget(curveTable);

    setWidget(dockWidgetContents);
}



