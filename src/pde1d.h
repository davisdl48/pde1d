#ifndef pde1d_H
#define pde1d_H

#include <QtGui/QMainWindow>
#include <QString>
#include <QtGui/QStandardItemModel>
#include "solvwidget.h"
#include "controls.h"
#include "ui_firstwindow.h"
#include "errtabdock.h"
#include "curvetabdock.h"
#include "curvesmodel.h"
#include <QList>

class pde1d : public QMainWindow, Ui_MainWindow
{
Q_OBJECT
public:
    pde1d();
    virtual ~pde1d();
    
    public slots:
       void setSize(int ivalue) ;
       void setCycles(double value=1.0);
       
       void reset();
       void setSpeed(double value = 1.0) {
	 if(value == speed) return;
	 speed = value;
	 dirty = true;
	 if( dx != 0.0) {
	 CFL = dt*speed/dx;
	 }
       }
       
       void setTimeStep(double value = 0.01) {
	 if(value == dt) return;
	 dt = value;
	 dirty = true;
	 if( dx != 0.0) {
	 CFL = dt*speed/dx;
	 }
       }
	
	 
       void setCFL(double value=1.0);
       void setEquation(int value=0);
       void setViscosity(double value=0.001);

       void setStop(int value) ;
       void setStep(int value) ;
       void setDelay(int value=10) ;
       void run();
       void stopRun();
       void saveImage();
       void savePlot();
       void addSolver(int index = 0);
       void removeSolver(int index = 0);
       //void removeSolver(QObject *obj = 0);
       void refresh() ;
	  	 
      
private:
  void replot(const char * title);
  void metrics();
  void updateNames();
  void addIt(SolvWidget * solver) ;
  
  int stop;
  int step;
  int delay;
  size_t N;
  double cycles;
  double dt;
  double speed;
  double CFL;
  int equation;
  double visc;
  QList<SolvWidget *> eeWidgets;
  CurveTabDock *curveTab;
  CurvesModel *curveModel;
  ErrTabDock *errTab;
  Controls *control; 
  QStandardItemModel * solvModel;
  int dockID;
  QwtLegend * legend;
};

#endif // pde1d_H
