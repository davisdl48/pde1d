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


#ifndef COLORDELEGATE_H
#define COLORDELEGATE_H

#include <qstyleditemdelegate.h>
#include <QComboBox>
#include <QList>
#include <QIcon>


class ColorDelegate : public QStyledItemDelegate
{
    Q_OBJECT
public:
    explicit ColorDelegate(QObject* parent = 0);
    virtual QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const;
    virtual void setEditorData(QWidget* editor, const QModelIndex& index) const;
    virtual void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const;
    virtual ~ColorDelegate();

    static const QIcon getIcon(QColor col) {
        QPixmap pmap(32,32);
        pmap.fill(col);
        return QIcon(pmap);
    }

    static QList<QColor> colors;
    static bool cinit;
    
    static bool initColors() {
        colors.append(Qt::red);
        colors.append(Qt::blue);
        colors.append(Qt::black);
        colors.append(Qt::darkCyan); /// #name
        colors.append(Qt::darkMagenta); /// #name
        colors.append(Qt::darkRed); /// #name
        colors.append(Qt::darkGreen); /// #name
        colors.append(Qt::darkBlue); /// #name
        colors.append(Qt::magenta);
        colors.append(Qt::lightGray); /// #name
        colors.append(Qt::cyan);
        colors.append(Qt::yellow);
        colors.append(Qt::white);
        colors.append(Qt::gray); /// darkGray
        colors.append(Qt::darkGray); /// gray #name
        colors.append(Qt::green); /// darkGreen
        colors.append(Qt::darkYellow); /// black #name
	return true;
    }

};

#endif // COLORDELEGATE_H
