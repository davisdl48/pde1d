/*
    To show a combobox editor in a table view - adapted from qt wiki
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


#include "itemdelegate.h"
#include <QComboBox>
#include <QStringList>
#include "myinputs.h"
#include "curvesmodel.h"
#include <iostream>
#include <qwt_symbol.h>


QList<NV> SymbolStyleDelegate::syms;
bool SymbolStyleDelegate::sinit = SymbolStyleDelegate::initSyms();

PenStyleDelegate::PenStyleDelegate(QObject* parent): QStyledItemDelegate(parent)
{
    
}

QWidget* PenStyleDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const
{
    QComboBox *cb = new QComboBox(parent);
    cb->addItem(tr("No Line"),Qt::NoPen);
    cb->addItem(tr("Solid Line"),Qt::SolidLine);
    cb->addItem(tr("Dashed"),Qt::DashLine);
    cb->addItem(tr("Dotted"),Qt::DotLine);
    cb->addItem(tr("Dash Dot"),Qt::DashDotLine);
    cb->addItem(tr("Dash Dot Dot"),Qt::DashDotDotLine);
    // cb->setCurrentIndex(index.data(Qt::EditRole));
    return cb;
}

void PenStyleDelegate::setEditorData(QWidget* editor, const QModelIndex& index) const
{
    int cbIndex;
    if(QComboBox *cb = qobject_cast<QComboBox *>(editor)) {
        // get the index of the text in the combobox that matches the current value of the item
        cbIndex = index.data(Qt::EditRole).toInt();
        // if it is valid, adjust the combobox
        if(cbIndex >= 0)  cb->setCurrentIndex(cbIndex);
    } else {
        QStyledItemDelegate::setEditorData(editor, index);
    }
}

void PenStyleDelegate::setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const
{
    if(QComboBox *cb = qobject_cast<QComboBox *>(editor)) {
        // save the current text of the combo box as the current value of the item
        model->setData(index, cb->currentIndex(), Qt::EditRole);
	//model->setData(index, cb->currentText(), Qt::DisplayRole);
    } else {
        QStyledItemDelegate::setModelData(editor, model, index);
    }
}


LineWidthDelegate::LineWidthDelegate(QObject* parent): QStyledItemDelegate(parent)
{
    
}

QWidget* LineWidthDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const
{
    MyDoubInput *dinp = new MyDoubInput(0.0,parent,0.0,100.0,0.1,2);
    return dinp;
}

void LineWidthDelegate::setEditorData(QWidget* editor, const QModelIndex& index) const
{
    bool ok;
    if(MyDoubInput *dinp = qobject_cast<MyDoubInput *>(editor)) {
        dinp->setValue(index.data(Qt::EditRole).toDouble(&ok));
    } else {
        QStyledItemDelegate::setEditorData(editor, index);
    }
}

void LineWidthDelegate::setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const
{
    MyDoubInput *dinp = qobject_cast<MyDoubInput *>(editor);
    if(dinp) {
      double lw = dinp->getValue();
      QVariant lwvar(lw);
      model->setData(index,lwvar, Qt::EditRole);
	
    } else {
        QStyledItemDelegate::setModelData(editor, model, index);
    }
}


SymbolStyleDelegate::SymbolStyleDelegate(QObject* parent): QStyledItemDelegate(parent)
{
    
}

QWidget* SymbolStyleDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const
{
    QComboBox *cb = new QComboBox(parent);
    QListIterator<NV> iter(syms);
    QwtSymbol symbol;
    QColor color;
    int size;
    QModelIndex ind;
    const CurvesModel *model = dynamic_cast<const CurvesModel*>(index.model());
    if( model == NULL) {
      symbol.setBrush(QBrush(Qt::white));
      symbol.setPen(QPen(Qt::black));
      symbol.setSize(QSize(10,10));
    }else{
      // get color and size from the model  
      color = model->getSymbolColor(index.row());
      size = model->getSymbolSize(index.row());
      symbol.setBrush(QBrush(Qt::white));
      symbol.setPen(QPen(color));
      symbol.setSize(QSize(size,size));
    }
    while(iter.hasNext()) {
      NV sym = iter.next();    
      //cb->addItem(getIcon(col),CurvesModel::colorName(col),variant);
      symbol.setStyle(sym.val);
      
      cb->addItem(getIcon(symbol),sym.name,sym.val);
    }
  
    //cb->addItem(tr("UserStyle"),15);
    return cb;
}

void SymbolStyleDelegate::setEditorData(QWidget* editor, const QModelIndex& index) const
{
    int cbIndex;
    if(QComboBox *cb = qobject_cast<QComboBox *>(editor)) {
        // get the index of the text in the combobox that matches the current value of the item
        cbIndex = cb->findData(index.data(Qt::EditRole).toInt());
        // if it is valid, adjust the combobox
        if(cbIndex >= -1)  {
	  cb->setCurrentIndex(cbIndex);
	}else{
	  cb->setCurrentIndex(0);
	}
	
    } else {
        QStyledItemDelegate::setEditorData(editor, index);
    }
}

void SymbolStyleDelegate::setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const
{ 
    int cbIndex;
    if(QComboBox *cb = qobject_cast<QComboBox *>(editor)) {
        cbIndex = cb->currentIndex();
        // save the current text of the combo box as the current value of the item
        model->setData(index, cb->itemData(cbIndex).toInt(), Qt::EditRole);
	//model->setData(index, cb->currentText(), Qt::DisplayRole);
    } else {
        QStyledItemDelegate::setModelData(editor, model, index);
    }
}


SymbolSizeDelegate::SymbolSizeDelegate(QObject* parent): QStyledItemDelegate(parent)
{
    
}

QWidget* SymbolSizeDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const
{
    MyIntInput *iinp = new MyIntInput(parent);
    return iinp;
}

void SymbolSizeDelegate::setEditorData(QWidget* editor, const QModelIndex& index) const
{

    if(MyIntInput *iinp = qobject_cast<MyIntInput *>(editor)) {
        iinp->setValue(index.data(Qt::EditRole).toInt());
    } else {
        QStyledItemDelegate::setEditorData(editor, index);
    }
}

void SymbolSizeDelegate::setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const
{
    MyIntInput *iinp = qobject_cast<MyIntInput *>(editor);
    if(iinp) {
      int ss = iinp->getValue();
      QVariant ssvar(ss);
      model->setData(index,ssvar, Qt::EditRole);
    } else {
        QStyledItemDelegate::setModelData(editor, model, index);
    }
}

