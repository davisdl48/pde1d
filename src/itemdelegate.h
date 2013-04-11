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


#ifndef ITEMDELEGATE_H
#define ITEMDELEGATE_H

#include <qstyleditemdelegate.h>
#include <QComboBox>
#include <QStringList>
#include <QVariant>
#include <qwt_symbol.h>

typedef struct { QString name; QwtSymbol::Style val; } NV;

class PenStyleDelegate : public QStyledItemDelegate
{
public:
    explicit PenStyleDelegate(QObject* parent = 0);
    virtual QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const;
    virtual void setEditorData(QWidget* editor, const QModelIndex& index) const;
    virtual void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const;
};

class LineWidthDelegate : public QStyledItemDelegate
{
public:
    explicit LineWidthDelegate(QObject* parent = 0);
    virtual QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const;
    virtual void setEditorData(QWidget* editor, const QModelIndex& index) const;
    virtual void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const;
};


class SymbolStyleDelegate : public QStyledItemDelegate
{
public:
    explicit SymbolStyleDelegate(QObject* parent = 0);
    virtual QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const;
    virtual void setEditorData(QWidget* editor, const QModelIndex& index) const;
    virtual void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const;
    
    static QList<NV> syms;
    static bool sinit;
    static bool initSyms() {
      NV pair;
      pair.name=tr("No Symbol");
      pair.val=QwtSymbol::NoSymbol;
        syms.append(pair);
      pair.name=tr("Ellipse");
      pair.val=QwtSymbol::Ellipse;
        syms.append (pair);
      pair.name=tr("Rectangle");
      pair.val=QwtSymbol::Rect;
        syms.append(pair);
      pair.name=tr("Diamond");
      pair.val=QwtSymbol::Diamond;
        syms.append(pair);
      pair.name=tr("Triangle");
      pair.val=QwtSymbol::Triangle;
        syms.append(pair);
      pair.name=tr("DTriangle");
      pair.val=QwtSymbol::DTriangle;
        syms.append(pair);
      pair.name=tr("UTriangle");
      pair.val=QwtSymbol::UTriangle;
        syms.append(pair);
      pair.name=tr("LTriangle");
      pair.val=QwtSymbol::LTriangle;
        syms.append(pair);
      pair.name=tr("RTriangle");
      pair.val=QwtSymbol::RTriangle;
        syms.append(pair);
      pair.name=tr("Cross");
      pair.val=QwtSymbol::Cross;
        syms.append(pair);
      pair.name=tr("XCross");
      pair.val=QwtSymbol::XCross;
        syms.append(pair);
      pair.name=tr("HLine");
      pair.val=QwtSymbol::HLine;
        syms.append(pair);
      pair.name=tr("VLine");
      pair.val=QwtSymbol::VLine;
        syms.append(pair);
      pair.name=tr("Star1");
      pair.val=QwtSymbol::Star1;
        syms.append(pair);
      pair.name=tr("Star2");
      pair.val=QwtSymbol::Star2;
        syms.append(pair);
      pair.name=tr("Hexagon");
      pair.val=QwtSymbol::Hexagon;
        syms.append(pair);
      pair.name= tr("User Style");
      pair.val=QwtSymbol::UserStyle;
      syms.append(pair);
         return true;
    }
    
    static QString symName(QwtSymbol::Style style) {
        QListIterator<NV> iter(syms);
    while(iter.hasNext()) {
      NV sym = iter.next();    
      //cb->addItem(getIcon(col),CurvesModel::colorName(col),variant);
      if( sym.val == style) return sym.name;
    }
    return "What??";
    }
      
};

class SymbolSizeDelegate : public QStyledItemDelegate
{
public:
    explicit SymbolSizeDelegate(QObject* parent = 0);
    virtual QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const;
    virtual void setEditorData(QWidget* editor, const QModelIndex& index) const;
    virtual void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const;
};

#endif

// ITEMDELEGATE_H
