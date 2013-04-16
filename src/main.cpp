#include <QtGui/QApplication>
#include "pde1d.h"


int main(int argc, char** argv)
{
  
  QApplication app(argc, argv);
  Q_INIT_RESOURCE(res); 
  QLabel l;
  l.setPixmap(QPixmap(":/images/linAdvect.png"));
  l.show();
  pde1d foo;
  foo.show();
  return app.exec();
  
  
}
