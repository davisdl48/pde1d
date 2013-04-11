/*
    Model for changing plot curve attributes
    Copyright (C) 2013  Doug Davis davisdl48@gmail.com

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


#ifndef CURVESMODEL_H
#define CURVESMODEL_H

#include <QAbstractItemModel>
#include "solvwidget.h"
#include <QPixmap>
#include <QPainter>


class CurvesModel : public QAbstractTableModel
{
Q_OBJECT
public:
    explicit CurvesModel(QObject* parent = 0);
    virtual QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const;
    virtual int columnCount(const QModelIndex& parent = QModelIndex()) const;
    virtual int rowCount(const QModelIndex& parent = QModelIndex()) const;
    virtual ~CurvesModel();
    virtual bool setData(const QModelIndex& index, const QVariant& value, int role = Qt::EditRole);
    virtual QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;
    virtual bool insertRows(int row, int count, const QModelIndex& parent = QModelIndex());
    virtual bool removeRows(int row, int count, const QModelIndex& parent = QModelIndex());
    virtual Qt::ItemFlags flags(const QModelIndex& index) const;

    void setSolvers(QList<SolvWidget *> *eeWidgets) ;

    void solverAdded(  ) ;

    void removeSolver( int row ) ;

    static QString colorName(QColor col) {
        if ( col == Qt::red ) return tr("Red");
        if ( col == Qt::blue ) return tr("Blue");
        if ( col == Qt::black ) return tr("Black");
        if ( col == Qt::green ) return tr("Green");
        if ( col == Qt::darkCyan ) return tr("DarkCyan");
        if ( col == Qt::darkMagenta ) return tr("DarkMagenta");
        if ( col == Qt::darkYellow ) return tr("DarkYellow");
        if ( col == Qt::darkGray ) return tr("DarkGray");
        if ( col == Qt::darkRed ) return tr("DarkRed");
        if ( col == Qt::darkGreen ) return tr("DarkGreen");
        if ( col == Qt::darkBlue ) return tr("DarkBlue");
        if ( col == Qt::magenta ) return tr("Magenta");
        if ( col == Qt::gray ) return tr("Gray");
        if ( col == Qt::lightGray ) return tr("LightGray");
        if ( col == Qt::cyan ) return tr("Cyan");
        if ( col == Qt::yellow ) return tr("Yellow");
        if ( col == Qt::white ) return tr("White");
        return col.name();
    }

    
    static QIcon penIcon(QPen pen) {
        QPixmap pmap(24,24);
	QColor col = pen.color();
	int bright = col.red() + 2*col.green() + col.blue();
        if(bright > 600 ) {
            pmap.fill(Qt::black);
        } else {
            pmap.fill(Qt::white);
        }
        QPainter painter(&pmap);
        painter.setPen(pen);
        QLineF line(0,11,23,11);
        painter.drawLine(line);
        return QIcon(pmap);
    }

signals:
  void newdata();
  
protected:
    QList<SolvWidget *> *solvers;
    QStringList penStyles;
};

#endif // CURVESMODEL_H
