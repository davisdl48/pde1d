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


#include "colordelegate.h"
#include "curvesmodel.h"
#include <QColor>
#include <iostream>
#include <QColorDialog>

QList<QColor> ColorDelegate::colors;
bool ColorDelegate::cinit = ColorDelegate::initColors();

ColorDelegate::ColorDelegate(QObject* parent): QStyledItemDelegate(parent)
{ 
}

QWidget* ColorDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const
{
  QComboBox *cb = new QComboBox(parent);
  QListIterator<QColor> iter(colors);
  while(iter.hasNext()) {
    QColor col = iter.next();
    QVariant variant = col;
    cb->addItem(getIcon(col),CurvesModel::colorName(col),variant);
  }
  cb->addItem(tr("Custom"));
  return cb;
}

void ColorDelegate::setEditorData(QWidget* editor, const QModelIndex& index) const
{
  if(QComboBox *cb = qobject_cast<QComboBox *>(editor)) { 
    int current = cb->findData(index.data(Qt::EditRole).value<QColor>(),Qt::EditRole );
    if(current > -1) {
      cb->setCurrentIndex(current);
    }else{
      cb->setCurrentIndex(index.row());
    }
    return;
  }
  QStyledItemDelegate::setEditorData(editor, index);
}

void ColorDelegate::setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const
{
  QColor col;
  if(QComboBox *cb = qobject_cast<QComboBox *>(editor)) {
     int current = cb->currentIndex();
     if(current == cb->count()-1 ) {
       col = QColorDialog::getColor(Qt::red,editor);
       if (col.isValid()) {  
	 model->setData(index,col,Qt::EditRole);
	 return;
       }
       return;
     }
     model->setData(index, cb->itemData(current,Qt::UserRole), Qt::EditRole);
  }
    QStyledItemDelegate::setModelData(editor, model, index);
}

ColorDelegate::~ColorDelegate()
{

}


