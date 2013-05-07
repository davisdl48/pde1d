/*
    QDockWidget with primary controls for pde1d
    Copyright (C) 2013  Davis Family davisdl48@gmail.com

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
    int gridrow;
    if (objectName().isEmpty())
        setObjectName(QString::fromUtf8("Controls"));
    //resize(304, 696);
    setAllowedAreas(Qt::LeftDockWidgetArea|Qt::RightDockWidgetArea);
    dockWidgetContents = new QWidget();
    dockWidgetContents->setObjectName(QString::fromUtf8("dockWidgetContents"));
    verticalLayout = new QVBoxLayout(dockWidgetContents);
    verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
    pdeLabel = new QLabel(tr("PDE"),dockWidgetContents);
    pdeLabel->setObjectName(QString::fromUtf8("pdeLabel"));

    verticalLayout->addWidget(pdeLabel);

    pdeBox = new QComboBox(dockWidgetContents);
    pdeBox->setObjectName(QString::fromUtf8("pdeBox"));

    verticalLayout->addWidget(pdeBox);

    //solveLabel = new QLabel(dockWidgetContents);
    //solveLabel->setObjectName(QString::fromUtf8("solveLabel"));

    //verticalLayout->addWidget(solveLabel);

    addSolvCombo = new QComboBox(dockWidgetContents);
    addSolvCombo->setObjectName(QString::fromUtf8("addSolvCombo"));

    verticalLayout->addWidget(addSolvCombo);

    

    cyclesLabel = new QLabel(dockWidgetContents);
    cyclesLabel->setObjectName(QString::fromUtf8("cyclesLabel"));

    verticalLayout->addWidget(cyclesLabel);

    cyclesInput = new MyDoubInput(dockWidgetContents);
    cyclesInput->setObjectName(QString::fromUtf8("cyclesInput"));
    cyclesInput->setMinimumWidth(30);
    cyclesInput->setMaximumWidth(222);
    cyclesInput->setValue(1);
    //cyclesInput->setMinimum(-9999);
    //cyclesInput->setDecimals(6);

    verticalLayout->addWidget(cyclesInput);

    gridrow = 0;
    gridLayout = new QGridLayout();
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
    
    speedLabel = new QLabel( tr("Speed"));
    speedLabel->setObjectName(tr("speedLabel"));
    gridLayout->addWidget(speedLabel, gridrow, 0, 1, Qt::AlignLeft);
    
    timeLabel = new QLabel( tr("Time Step"));
    timeLabel->setObjectName(tr("timeLabel"));
    gridLayout->addWidget(timeLabel, gridrow++, 1, 1, Qt::AlignLeft);
    
    speedInput = new MyDoubInput(dockWidgetContents);
    speedInput->setObjectName(tr("speedInput"));
    speedInput->setMinimumWidth(30);
    speedInput->setMaximumWidth(111);
    gridLayout->addWidget(speedInput, gridrow, 0, 1, Qt::AlignLeft);
    
    timeInput = new MyDoubInput(dockWidgetContents);
    timeInput->setObjectName(tr("timeInput"));
    timeInput->setMinimumWidth(30);
    timeInput->setMaximumWidth(111);
    gridLayout->addWidget(timeInput, gridrow++, 1, 1, Qt::AlignLeft);
    
    
    cflLabel = new QLabel(dockWidgetContents);
    cflLabel->setObjectName(QString::fromUtf8("cflLabel"));
    gridLayout->addWidget(cflLabel, gridrow, 0, 1, Qt::AlignLeft);

    viscLabel = new QLabel(dockWidgetContents);
    viscLabel->setObjectName(QString::fromUtf8("viscLabel"));
    gridLayout->addWidget(viscLabel, gridrow++, 1, 1, Qt::AlignLeft);

    cflInput = new MyDoubInput(dockWidgetContents);
    cflInput->setObjectName(QString::fromUtf8("cflInput"));
    cflInput->setMinimumWidth(30);
    cflInput->setMaximumWidth(111);
    cflInput->setValue(1);
    //cflInput->setMinimum(-9999);
    //cflInput->setDecimals(4);

    gridLayout->addWidget(cflInput, gridrow, 0, 1, Qt::AlignLeft);

    viscInput = new MyDoubInput(dockWidgetContents);
    viscInput->setObjectName(QString::fromUtf8("viscInput"));
    viscInput->setMinimumWidth(30);
    viscInput->setMaximumWidth(111);
    //viscInput->setMinimum(-1);
    //viscInput->setSingleStep(0.0001);
    //viscInput->setDecimals(4);

    gridLayout->addWidget(viscInput, gridrow++, 1, 1, Qt::AlignLeft);

    sizeLabel = new QLabel(dockWidgetContents);
    sizeLabel->setObjectName(QString::fromUtf8("sizeLabel"));
    gridLayout->addWidget(sizeLabel, gridrow++, 0, 1, 2, Qt::AlignLeft);

    intNumberOfPoints = new MyIntInput(dockWidgetContents);
    intNumberOfPoints->setObjectName(QString::fromUtf8("intNumberOfPoints"));
    intNumberOfPoints->setMinimumWidth(30);
    intNumberOfPoints->setMaximumWidth(111);
    intNumberOfPoints->setValue(100);
    //intNumberOfPoints->setMinimum(4);   
    gridLayout->addWidget(intNumberOfPoints, gridrow++, 0, 1, 2, Qt::AlignLeft);

    stepsLabel = new QLabel(dockWidgetContents);
    stepsLabel->setObjectName(QString::fromUtf8("stepsLabel"));
    gridLayout->addWidget(stepsLabel, gridrow, 0, 1, Qt::AlignLeft);

    incrementLabel = new QLabel(dockWidgetContents);
    incrementLabel->setObjectName(QString::fromUtf8("incrementLabel"));
    gridLayout->addWidget(incrementLabel, gridrow++, 1, 1, Qt::AlignLeft);

    intTimeSteps = new MyIntInput(dockWidgetContents);
    intTimeSteps->setObjectName(QString::fromUtf8("intTimeSteps"));
    intTimeSteps->setMinimumWidth(30);
    intTimeSteps->setMaximumWidth(111);
    intTimeSteps->setValue(100);
    //intTimeSteps->setMinimum(0);

    gridLayout->addWidget(intTimeSteps, gridrow, 0, 1, Qt::AlignLeft);

    intPlotIncrement = new MyIntInput(dockWidgetContents);
    intPlotIncrement->setObjectName(QString::fromUtf8("intPlotIncrement"));
    intPlotIncrement->setMinimumWidth(30);
    intPlotIncrement->setMaximumWidth(111);
    intPlotIncrement->setValue(1);
    //intPlotIncrement->setMinimum(1);

    gridLayout->addWidget(intPlotIncrement, gridrow++, 1, 1, Qt::AlignLeft);

    label = new QLabel(dockWidgetContents);
    label->setObjectName(QString::fromUtf8("label"));

    gridLayout->addWidget(label, gridrow, 0, 1, Qt::AlignLeft);

    plotDelayInput = new MyIntInput(dockWidgetContents);
    plotDelayInput->setObjectName(QString::fromUtf8("plotDelayInput"));
    plotDelayInput->setMinimumWidth(30);
    plotDelayInput->setMaximumWidth(111);

    gridLayout->addWidget(plotDelayInput, gridrow++, 1, 1, Qt::AlignLeft);


    verticalLayout->addLayout(gridLayout);
    
    horizontalLayout2 = new QHBoxLayout();
    horizontalLayout2->setObjectName(QString::fromUtf8("horizontalLayout2"));
    
    saveImageButton = new QPushButton(dockWidgetContents);
    saveImageButton->setObjectName(QString::fromUtf8("saveImageButton"));
    horizontalLayout2->addWidget(saveImageButton);
    
    savePlotButton = new QPushButton(dockWidgetContents);
    savePlotButton->setObjectName(QString::fromUtf8("savePlotButton"));
    horizontalLayout2->addWidget(savePlotButton);

    verticalLayout->addLayout(horizontalLayout2);

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
    //solveLabel->setBuddy(addSolvCombo);
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
    pdeLabel->setText(QApplication::translate("Controls", "PDE", 0, QApplication::UnicodeUTF8));
    //solveLabel->setText(QApplication::translate("Controls", "Add Solver", 0, QApplication::UnicodeUTF8));
    sizeLabel->setText(QApplication::translate("Controls", "Number of Points", 0, QApplication::UnicodeUTF8));
    //intNumberOfPoints->setLabel(QString());
    cyclesLabel->setText(QApplication::translate("Controls", "Number of Cycles", 0, QApplication::UnicodeUTF8));
    cflLabel->setText(QApplication::translate("Controls", "CFL", 0, QApplication::UnicodeUTF8));
    viscLabel->setText(QApplication::translate("Controls", "Viscosity", 0, QApplication::UnicodeUTF8));
    //cflInput->setLabel(QString());
    stepsLabel->setText(QApplication::translate("Controls", "Time Steps", 0, QApplication::UnicodeUTF8));
    incrementLabel->setText(QApplication::translate("Controls", "Plot Incr.", 0, QApplication::UnicodeUTF8));
    label->setText(QApplication::translate("Controls", "Plot Delay", 0, QApplication::UnicodeUTF8));
    saveImageButton->setText(QApplication::translate("Controls", "Save Image", 0, QApplication::UnicodeUTF8));
    savePlotButton->setText(QApplication::translate("Controls", "Save Plot", 0, QApplication::UnicodeUTF8));
    animationCheck->setText(QApplication::translate("Controls", "Save Anim. Frames", 0, QApplication::UnicodeUTF8));
    checkBox->setText(QApplication::translate("Controls", "Save Error Data", 0, QApplication::UnicodeUTF8));
    resetButton->setText(QApplication::translate("Controls", "Reset", 0, QApplication::UnicodeUTF8));
    runButton->setText(QApplication::translate("Controls", "Run", 0, QApplication::UnicodeUTF8));
    stopButton->setText(QApplication::translate("Controls", "Stop", 0, QApplication::UnicodeUTF8));
}

