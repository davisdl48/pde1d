/*
    Model for changing plot curve attributes
    Copyright (C) 2013  Doug Davis <davisdl48@gmail.com>

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

#include "curvesmodel.h"
#include "solvwidget.h"
#include "colordelegate.h"
#include "itemdelegate.h"
#include <iostream>
#include <qwt_symbol.h>

CurvesModel::CurvesModel(QObject* parent): QAbstractTableModel(parent)
{
    penStyles << "No Line"<<"Solid Line"<<"Dashed"<<"Dotted"<<"Dash Dot"<<"Dash Dot Dot";
}

QVariant CurvesModel::data(const QModelIndex& index, int role) const
{
    SolvWidget* solv =  solvers->at(index.row());
    QwtPlotCurve * curve = solv->curve;
    QColor col = curve->pen().color();
    int bright = col.red() + 2*col.green() + col.blue();
    QwtSymbol::Style symSty;
    switch(index.column()) {
    case 0:// title
        switch(role) {
        case Qt::DisplayRole :
        case Qt::EditRole:
            return solv->getTitle();
        }
        break;
    case 1:// line color
        switch(role) {
        case Qt::BackgroundRole:
            return QBrush(col);
        case Qt::ForegroundRole:
            if(bright > 600 ) return QBrush(Qt::darkGray);
            return QBrush(Qt::white);
        case Qt::DisplayRole:
            return QVariant(colorName(col));
        case Qt::EditRole:
        case Qt::UserRole:
            return QVariant(col);
        }
        break;
    case 2:// line style
        switch(role) {
        case Qt::DecorationRole:
            return penIcon(solv->curve->pen());
        case Qt::DisplayRole :
            return penStyles[solv->curve->pen().style()];
        case Qt::EditRole :
            return solv->curve->pen().style();
        }
        break;
    case 3:// line width
        switch(role) {
        case Qt::DecorationRole:
            return penIcon(solv->curve->pen());
        case Qt::DisplayRole :
            return QString("%1").arg(solv->curve->pen().widthF(),4,'f',2);
        case Qt::EditRole :
            return solv->curve->pen().widthF();
        }
        break;
    case 4:// symbol style
        if(curve->symbol() == NULL) {
            symSty = QwtSymbol::NoSymbol;
        } else {
            symSty = curve->symbol()->style();
        }
        switch(role) {
        case Qt::EditRole :
            return QVariant(symSty);
        case Qt::DisplayRole :
            return SymbolStyleDelegate::symName(symSty);
        }
        break;
    case 5:// symbol color
        if(curve->symbol() != NULL && curve->symbol()->style() != QwtSymbol::NoSymbol) {
            col= curve->symbol()->pen().color();
        } // default to curve pen color
        switch(role) {
        case Qt::BackgroundRole:
            return QBrush(col);
        case Qt::ForegroundRole:
            if(bright > 600 ) return QBrush(Qt::darkGray);
            return QBrush(Qt::white);
        case Qt::DisplayRole :
            if(curve->symbol() == NULL || curve->symbol()->style() == QwtSymbol::NoSymbol) {
                return ("No Symbol");
            }
            return QVariant(colorName(col)) ;
        case Qt::EditRole:
        case Qt::UserRole:
            return QVariant(col);
        }
        break;
    case 6: // symbol size
        switch(role) {
            //case Qt::DecorationRole:
            //return penIcon(solv->curve->pen());
        case Qt::DisplayRole :
	  if(curve->symbol() != NULL && curve->symbol()->style() != QwtSymbol::NoSymbol ) {
            return QString("%1").arg(curve->symbol()->size().width());
	  }
	  break;
        case Qt::EditRole :
        case Qt::UserRole:
	  if(curve->symbol() != NULL && curve->symbol()->style() != QwtSymbol::NoSymbol ) {
            return curve->symbol()->size().width();
	  }
	  return 7;
        }
    }
    return QVariant();
}

int CurvesModel::columnCount(const QModelIndex& parent) const
{
    return 7;
}

int CurvesModel::rowCount(const QModelIndex& parent) const
{
    if(solvers->empty()) return 0;
    return solvers->size();
}

CurvesModel::~CurvesModel()
{

}

//Q_DECLARE_METATYPE (Qt::PenStyle)

bool CurvesModel::setData(const QModelIndex& index, const QVariant& value, int role)
{
    SolvWidget* solv =  solvers->at(index.row());
    QwtPlotCurve* curve = solv->curve;
    QPen pen = curve->pen(); /// This makes a copy of Pen
    QColor col;
    Qt::PenStyle ps;
    QwtSymbol *symbol;
    QwtSymbol::Style symbolStyle;
    QBrush brush;
    QSize symSize;

    bool ok;
    double lw;
    int ss;
    if( role != Qt::EditRole ) return false;
    switch(index.column()) {
    case 0:// title
        if(value.toString().size() > 1) {
            QString newTitle = value.toString();
            if( newTitle != solv->getTitle() ) {
                solv->setTitle(value.toString());
                emit newdata();
                return true;
            }
        }
        return false;
    case 1:// line color
        col = value.value<QColor>();
        if(pen.color() != col ) {
            pen.setColor(col);
            curve->setPen(pen);
            solv->setColor(col);
            emit newdata();
            return true;
        }
        return false;
    case 2:// line style
        ps = Qt::PenStyle(value.toInt());
        if( ps != pen.style()) {
            pen.setStyle(ps);
            curve->setPen(pen);
            emit newdata();
            return true;
        }
        return false;
    case 3:// line width
        lw = value.toDouble(&ok);
        if(ok && (lw != curve->pen().widthF())) {
            pen.setWidthF(lw);
            curve->setPen(pen);
            emit newdata();
            return true;
        }
        return false;
    case 4:// symbol style
        symbolStyle = QwtSymbol::Style(value.toInt());
        if(curve->symbol() == NULL) {
            if(symbolStyle == QwtSymbol::NoSymbol ) return false;
	    pen.setStyle(Qt::SolidLine);
            symbol = new QwtSymbol(symbolStyle,QBrush(),pen,QSize(7,7));//Qt::transparent
        } else {
            if(curve->symbol()->style() == symbolStyle) return false;
            if(symbolStyle == QwtSymbol::NoSymbol) {
                symbol = NULL;
            } else {
                symbol = new QwtSymbol(symbolStyle,curve->symbol()->brush(),curve->symbol()->pen(),curve->symbol()->size());
            }
        }
        curve->setSymbol(symbol);
        emit newdata();
        return true;
    case 5:// symbol color
        col = value.value<QColor>();
        if(curve->symbol() == NULL || curve->symbol()->style() == QwtSymbol::NoSymbol ) return false;
        if(curve->symbol()->pen().color() == col ) return false;
        pen = curve->symbol()->pen();
        pen.setColor(col);
        symbol = new QwtSymbol(curve->symbol()->style(),curve->symbol()->brush(),pen,curve->symbol()->size());
        curve->setSymbol(symbol);
        emit newdata();
        return true;
    case 6: // symbol size
        ss = value.toInt();
	if(ss < 1) return false;
        if(curve->symbol() == NULL || curve->symbol()->style() == QwtSymbol::NoSymbol ) return false;
        if(curve->symbol()->size().width() == ss ) return false;
        symbol = new QwtSymbol(curve->symbol()->style(),curve->symbol()->brush(),curve->symbol()->pen(),QSize(ss,ss));
        curve->setSymbol(symbol);
        emit newdata();
        return false;
    }
    return QAbstractItemModel::setData(index, value, role);
}

QVariant CurvesModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role != Qt::DisplayRole) return QVariant();
    if(orientation == Qt::Vertical) {
        return QVariant(solvers->at(section)->getTitle());
    }
    switch(section) {
    case 0:
        return  tr("Solver");
    case 1:
        return  tr("Color");
    case 2:
        return tr("Line Style");
    case 3:
        return tr("Line Width");
    case 4:
        return tr("Symbol");
    case 5:
        return tr("Sym Color");
    case 6:
        return tr("Sym Size");
    }
    return QAbstractItemModel::headerData(section, orientation, role);
}

bool CurvesModel::insertRows(int row, int count, const QModelIndex& parent)
{
    return QAbstractItemModel::insertRows(row, count, parent);
}

bool CurvesModel::removeRows(int row, int count, const QModelIndex& parent)
{
    return QAbstractItemModel::removeRows(row, count, parent);
}

Qt::ItemFlags CurvesModel::flags(const QModelIndex& index) const
{
    SolvWidget* solv =  solvers->at(index.row());
    QwtPlotCurve* curve = solv->curve;
    if(index.column() < 5) return Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsSelectable;
    if(index.column() == 5 || index.column() == 6) {
        if(curve->symbol() != NULL && curve->symbol()->style() != QwtSymbol::NoSymbol) {
            return Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsSelectable;
        } else {
            return Qt::NoItemFlags;
        }
    }
    return QAbstractItemModel::flags(index);
}

void CurvesModel::removeSolver(int row) {
}

void CurvesModel::solverAdded() {
    int nlast=solvers->size()-1;
    beginInsertRows(QModelIndex(),nlast,nlast);
    endInsertRows();
}

void CurvesModel::setSolvers(QList< SolvWidget* >  *eeWidgets) {
    beginResetModel();
    solvers = eeWidgets;
    endResetModel();
}

QColor CurvesModel::getSymbolColor(int row) const { 
    SolvWidget* solv =  solvers->at(row);
    QwtPlotCurve* curve = solv->curve;
    if(curve->symbol() != NULL && curve->symbol()->style() != QwtSymbol::NoSymbol) {
            return curve->symbol()->pen().color();
    }
    return curve->pen().color();
}

int CurvesModel::getSymbolSize(int row) const {
    SolvWidget* solv =  solvers->at(row);
    QwtPlotCurve* curve = solv->curve;
    if(curve->symbol() != NULL && curve->symbol()->style() != QwtSymbol::NoSymbol) {
            return curve->symbol()->size().width();
    }
    return 7;
}


