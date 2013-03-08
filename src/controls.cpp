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


#include "controls.h"

Controls::Controls ( const QString& title, QWidget* parent, Qt::WindowFlags flags ) : QDockWidget ( parent, flags )
{
  setWindowTitle(tr("Controls"));
  setupUi ();
}


Controls::Controls ( const Controls& other )
{

}

Controls::~Controls()
{

}

Controls& Controls::operator= ( const Controls& other )
{
    return *this;
}

bool Controls::operator== ( const Controls& other ) const
{
return this == &other;
}
void Controls::setupUi() {
    if (objectName().isEmpty())
        setObjectName(QString::fromUtf8("Controls"));
    //resize(304, 696);
    setAllowedAreas(Qt::LeftDockWidgetArea|Qt::RightDockWidgetArea);
    dockWidgetContents = new QWidget();
    dockWidgetContents->setObjectName(QString::fromUtf8("dockWidgetContents"));
    verticalLayout = new QVBoxLayout(dockWidgetContents);
    verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
    pdeLable = new QLabel(tr("PDE"),dockWidgetContents);
    pdeLable->setObjectName(QString::fromUtf8("pdeLable"));

    verticalLayout->addWidget(pdeLable);

    pdeBox = new QComboBox(dockWidgetContents);
    pdeBox->setObjectName(QString::fromUtf8("pdeBox"));

    verticalLayout->addWidget(pdeBox);

    solveLabel = new QLabel(dockWidgetContents);
    solveLabel->setObjectName(QString::fromUtf8("solveLabel"));

    verticalLayout->addWidget(solveLabel);

    addSolvCombo = new QComboBox(dockWidgetContents);
    addSolvCombo->setObjectName(QString::fromUtf8("addSolvCombo"));

    verticalLayout->addWidget(addSolvCombo);

    sizeLabel = new QLabel(dockWidgetContents);
    sizeLabel->setObjectName(QString::fromUtf8("sizeLabel"));

    verticalLayout->addWidget(sizeLabel);

    intNumberOfPoints = new MyIntInput(dockWidgetContents);
    intNumberOfPoints->setObjectName(QString::fromUtf8("intNumberOfPoints"));
    intNumberOfPoints->setValue(100);
    //intNumberOfPoints->setMinimum(4);

    verticalLayout->addWidget(intNumberOfPoints);

    cyclesLabel = new QLabel(dockWidgetContents);
    cyclesLabel->setObjectName(QString::fromUtf8("cyclesLabel"));

    verticalLayout->addWidget(cyclesLabel);

    cyclesInput = new MyDoubInput(dockWidgetContents);
    cyclesInput->setObjectName(QString::fromUtf8("cyclesInput"));
    cyclesInput->setValue(1);
    //cyclesInput->setMinimum(-9999);
    //cyclesInput->setDecimals(6);

    verticalLayout->addWidget(cyclesInput);

    gridLayout = new QGridLayout();
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
    cflLabel = new QLabel(dockWidgetContents);
    cflLabel->setObjectName(QString::fromUtf8("cflLabel"));

    gridLayout->addWidget(cflLabel, 0, 0, 1, 1);

    viscLabel = new QLabel(dockWidgetContents);
    viscLabel->setObjectName(QString::fromUtf8("viscLabel"));

    gridLayout->addWidget(viscLabel, 0, 1, 1, 1);

    cflInput = new MyDoubInput(dockWidgetContents);
    cflInput->setObjectName(QString::fromUtf8("cflInput"));
    cflInput->setValue(1);
    //cflInput->setMinimum(-9999);
    //cflInput->setDecimals(4);

    gridLayout->addWidget(cflInput, 1, 0, 1, 1);

    viscInput = new MyDoubInput(dockWidgetContents);
    viscInput->setObjectName(QString::fromUtf8("viscInput"));
    //viscInput->setMinimum(-1);
    //viscInput->setSingleStep(0.0001);
    //viscInput->setDecimals(4);

    gridLayout->addWidget(viscInput, 1, 1, 1, 1);

    stepsLabel = new QLabel(dockWidgetContents);
    stepsLabel->setObjectName(QString::fromUtf8("stepsLabel"));

    gridLayout->addWidget(stepsLabel, 2, 0, 1, 1);

    incrementLabel = new QLabel(dockWidgetContents);
    incrementLabel->setObjectName(QString::fromUtf8("incrementLabel"));

    gridLayout->addWidget(incrementLabel, 2, 1, 1, 1);

    intTimeSteps = new MyIntInput(dockWidgetContents);
    intTimeSteps->setObjectName(QString::fromUtf8("intTimeSteps"));
    intTimeSteps->setValue(100);
    //intTimeSteps->setMinimum(0);

    gridLayout->addWidget(intTimeSteps, 3, 0, 1, 1);

    intPlotIncrement = new MyIntInput(dockWidgetContents);
    intPlotIncrement->setObjectName(QString::fromUtf8("intPlotIncrement"));
    intPlotIncrement->setValue(1);
    //intPlotIncrement->setMinimum(1);

    gridLayout->addWidget(intPlotIncrement, 3, 1, 1, 1);

    label = new QLabel(dockWidgetContents);
    label->setObjectName(QString::fromUtf8("label"));

    gridLayout->addWidget(label, 4, 0, 1, 1);

    plotDelayInput = new MyIntInput(dockWidgetContents);
    plotDelayInput->setObjectName(QString::fromUtf8("plotDelayInput"));

    gridLayout->addWidget(plotDelayInput, 4, 1, 1, 1);


    verticalLayout->addLayout(gridLayout);

    savePlotButton = new QPushButton(dockWidgetContents);
    savePlotButton->setObjectName(QString::fromUtf8("savePlotButton"));

    verticalLayout->addWidget(savePlotButton);

    animationCheck = new QCheckBox(dockWidgetContents);
    animationCheck->setObjectName(QString::fromUtf8("animationCheck"));

    verticalLayout->addWidget(animationCheck);

    checkBox = new QCheckBox(dockWidgetContents);
    checkBox->setObjectName(QString::fromUtf8("checkBox"));

    verticalLayout->addWidget(checkBox);

    verticalSpacer_2 = new QSpacerItem(237, 125, QSizePolicy::Minimum, QSizePolicy::Expanding);

    verticalLayout->addItem(verticalSpacer_2);

    horizontalLayout = new QHBoxLayout();
    horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
    resetButton = new QPushButton(dockWidgetContents);
    resetButton->setObjectName(QString::fromUtf8("resetButton"));

    horizontalLayout->addWidget(resetButton);

    runButton = new QPushButton(dockWidgetContents);
    runButton->setObjectName(QString::fromUtf8("runButton"));

    horizontalLayout->addWidget(runButton);

    stopButton = new QPushButton(dockWidgetContents);
    stopButton->setObjectName(QString::fromUtf8("stopButton"));

    horizontalLayout->addWidget(stopButton);


    verticalLayout->addLayout(horizontalLayout);

    setWidget(dockWidgetContents);
#ifndef QT_NO_SHORTCUT
    solveLabel->setBuddy(addSolvCombo);
    sizeLabel->setBuddy(intNumberOfPoints);
    cyclesLabel->setBuddy(cyclesInput);
    cflLabel->setBuddy(cflInput);
    stepsLabel->setBuddy(intTimeSteps);
    incrementLabel->setBuddy(intPlotIncrement);
    label->setBuddy(plotDelayInput);
#endif // QT_NO_SHORTCUT

    retranslateUi();

    //QMetaObject::connectSlotsByName(Controls);
}


void Controls::retranslateUi() {
    setWindowTitle(QApplication::translate("Controls", "Controls", 0, QApplication::UnicodeUTF8));
    pdeLable->setText(QApplication::translate("Controls", "PDE", 0, QApplication::UnicodeUTF8));
    solveLabel->setText(QApplication::translate("Controls", "Add Solver", 0, QApplication::UnicodeUTF8));
    sizeLabel->setText(QApplication::translate("Controls", "Number of Points", 0, QApplication::UnicodeUTF8));
    //intNumberOfPoints->setLabel(QString());
    cyclesLabel->setText(QApplication::translate("Controls", "Number of Cycles", 0, QApplication::UnicodeUTF8));
    cflLabel->setText(QApplication::translate("Controls", "CFL", 0, QApplication::UnicodeUTF8));
    viscLabel->setText(QApplication::translate("Controls", "Viscosity", 0, QApplication::UnicodeUTF8));
    //cflInput->setLabel(QString());
    stepsLabel->setText(QApplication::translate("Controls", "Time Steps", 0, QApplication::UnicodeUTF8));
    incrementLabel->setText(QApplication::translate("Controls", "Plot Incr.", 0, QApplication::UnicodeUTF8));
    label->setText(QApplication::translate("Controls", "Plot Delay", 0, QApplication::UnicodeUTF8));
    savePlotButton->setText(QApplication::translate("Controls", "Save Image", 0, QApplication::UnicodeUTF8));
    animationCheck->setText(QApplication::translate("Controls", "Save Anim. Frames", 0, QApplication::UnicodeUTF8));
    checkBox->setText(QApplication::translate("Controls", "Save Error Data", 0, QApplication::UnicodeUTF8));
    resetButton->setText(QApplication::translate("Controls", "Reset", 0, QApplication::UnicodeUTF8));
    runButton->setText(QApplication::translate("Controls", "Run", 0, QApplication::UnicodeUTF8));
    stopButton->setText(QApplication::translate("Controls", "Stop", 0, QApplication::UnicodeUTF8));
}

