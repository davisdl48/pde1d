#include <QtGui/QApplication>
#include "pde1d.h"


int main(int argc, char** argv)
{
  
  QApplication app(argc, argv);
  pde1d foo;
  foo.show();
  return app.exec();
  
  
}
