#ifndef pde1d_H
#define pde1d_H

#include <QtGui/QMainWindow>
#include <QString>
#include <QtGui/QStandardItemModel>
#include "solvwidget.h"
#include "controls.h"
#include "ui_firstwindow.h"
#include "errtabdock.h"
#include <vector>

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
  double CFL;
  int equation;
  double visc;
  std::vector<SolvWidget *> eeWidgets;
  ErrTabDock *errTab;
  Controls *control; 
  QStandardItemModel * solvModel;
  int dockID;
};

#endif // pde1d_H
