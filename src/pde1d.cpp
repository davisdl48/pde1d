#include "pde1d.h"
#include "solvwidget.h"
#include "eewidget.h"
#include "leastsqrwidget.h"
#include "simpimpwidget.h"
#include "femwidget.h"
#include "impwidget.h"
#include "rkwidget.h"

#include <QtGui/QLabel>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QAction>
#include <QtGui/QPen>
#include <QTimer>
#include <QFileDialog>
#include <QRegExp>
#include <QtGui/QImage>
#include <qwt_plot.h>
#include <qwt_legend.h>
#include <qwt_plot_curve.h>
#include <iostream>
#include <sstream>
#include <iomanip>


pde1d::pde1d() : QMainWindow(), Ui_MainWindow()
{
    setupUi ( this );

    N = 100;
    cycles = 1.0;
    CFL = 1.0;
    stop = 100;
    step = 1;
    dockID = 0;

    addSolver ( 0 );
    setCorner(Qt::BottomLeftCorner,Qt::LeftDockWidgetArea);

    qwtPlot->setAxisScale ( 0, -1.2, 1.2 );
    qwtPlot->setAxisScale ( 2, 0, 6.2831853 );
    qwtPlot->setCanvasBackground ( Qt::white );

    errTab = new ErrTabDock ( tr ( "Metrics" ), this );
    errTab->setFeatures ( QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetMovable );
    addDockWidget ( Qt::BottomDockWidgetArea, errTab );


    control = new Controls ( tr ( "Controls" ), this );
    //qwtPlot->setBackgroundRole(QPalette::Base);
    //alphaLineEdit->setValidator(&qv1);
    connect ( control->addSolvCombo, SIGNAL ( activated ( int ) ), this, SLOT ( addSolver ( int ) ) );
    //connect ( control->removeSolvCombo, SIGNAL ( activated ( int ) ), this, SLOT ( removeSolver ( int ) ) );
    connect ( control->intNumberOfPoints, SIGNAL ( valueChanged ( int ) ), this, SLOT ( setSize ( int ) ) );
    //cyclesInput->setValue(1.0);
    connect ( control->cyclesInput, SIGNAL ( valueChanged ( double ) ), this, SLOT ( setCycles ( double ) ) );
    //cflInput->setValue(1.0);
    connect ( control->cflInput, SIGNAL ( valueChanged ( double ) ), this, SLOT ( setCFL ( double ) ) );
    //connect(alphaLineEdit,SIGNAL(editingFinished()),this,SLOT(setAlpha()));
    connect ( control->intTimeSteps, SIGNAL ( valueChanged ( int ) ), this, SLOT ( setStop ( int ) ) );
    connect ( control->intPlotIncrement, SIGNAL ( valueChanged ( int ) ), this, SLOT ( setStep ( int ) ) );
    connect ( control->savePlotButton, SIGNAL ( clicked ( bool ) ), this, SLOT ( saveImage() ) );
    //connect ( control->plotDelayInput,SIGNAL( valueChanged ( int ) , this, SLOT( ..... ));
    connect ( control->resetButton, SIGNAL ( clicked ( bool ) ), this, SLOT ( reset() ) );
    connect ( control->runButton, SIGNAL ( clicked ( bool ) ), this, SLOT ( run() ) );
    connect ( control->stopButton, SIGNAL ( clicked ( bool ) ), this, SLOT ( stopRun() ) );
    setTabPosition(Qt::LeftDockWidgetArea |Qt::RightDockWidgetArea, QTabWidget::West);
    control->setFeatures ( QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetMovable );
    addDockWidget ( Qt::LeftDockWidgetArea, control );
    replot ( "Initial Condition" );
    savePlot ( "Initial.png" );


    //lw->setup();
    //void setup(const size_t size=100,const double cycles=1.0,const double waveSpeed =1.0,const double cfl =1.0,const double alpha=0.5);
    // size, cycles,  CFL

}

pde1d::~pde1d()
{}

void pde1d::replot ( const char * title )
{
    if ( eeWidgets.empty() ) return;

    qwtPlot->detachItems();
    qwtPlot->insertLegend ( new QwtLegend(), QwtPlot::BottomLegend );

    QwtPlotCurve *c1;

    c1 = new QwtPlotCurve ( "Ideal Values" );
    c1->setSamples ( eeWidgets[0]->getX(), eeWidgets[0]->getIdeal(), N );
    c1->attach ( qwtPlot );

    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        c1 = new QwtPlotCurve ( eeWidgets[i]->getTitle() );
        c1->setSamples ( eeWidgets[i]->getX(), eeWidgets[i]->getU(), N );
        c1->setPen ( QPen ( QBrush ( eeWidgets[i]->getColor() ), 1, Qt::DashLine ) );
        //c1->pen().setStyle(Qt::DashLine);
        c1->attach ( qwtPlot );
    }

    qwtPlot->setTitle ( title );
    qwtPlot->replot();
}

void pde1d::metrics()
{
    double *ideal;
    double *sim;
    double maxerr, rmserr, maxval, minval, totvar;
    double u, oldu, du;
    double err;
    QStringList colname;
    QTableWidgetItem oldname;
    QString ti;
    SolvWidget *solv;
    if ( eeWidgets.empty() ) return;
    if ( ( size_t ) errTab->errTabWid->columnCount() != eeWidgets.size() ) errTab->errTabWid->setColumnCount ( eeWidgets.size() );
    if ( errTab->errTabWid->rowCount() < 5 ) errTab->errTabWid->setRowCount ( 5 );
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        solv = eeWidgets[i];
        //std::cout << "pde1d::metrics solver " << solv->getTitle().toLocal8Bit().data() << "  " << i << "  " << N << std::endl;
        ideal = solv->getIdeal();
        sim = solv->getU();
        maxerr = 0.0;
        rmserr = 0.0;
        maxval = -1e32;
        minval = 1e32;
        totvar = 0.0;
        oldu = sim[0];
        for ( size_t j = 0; j < N; j++ ) {
            u = sim[j];
            maxval = ( u > maxval ) ? u : maxval;
            minval = ( u < minval ) ? u : minval;
            du = u - oldu;
            du = ( du > 0.0 ) ? du : -du;
            totvar += du;
            oldu = u;
            err = ideal[j] - u;
            rmserr += err * err;
            err = ( err > 0.0 ) ? err : -err;
            maxerr = ( err > maxerr ) ? err : maxerr;
        }
        du = sim[0] - oldu;
        du = ( du > 0.0 ) ? du : -du;
        totvar += du;
        rmserr = sqrt ( rmserr / N );
        //std::cout << maxerr << '\t' << rmserr << '\t' << maxval << '\t' << minval << '\t' << totvar << std::endl;
        colname.append ( eeWidgets[i]->getTitle() );
        // add data to errTab->errTab
        //maxerr
        //std::cout << "maxerr =  " << maxerr << std::endl;
        ti = QString ( "%1" ).arg ( maxerr );
        if ( errTab->errTabWid->item ( 0, i ) == 0 ) {
            errTab->errTabWid->setItem ( 0, i, new QTableWidgetItem ( ti ) );
        } else {
            errTab->errTabWid->item ( 0, i )->setText ( ti );
        }
        //rmserr
        //std::cout << "rmserr =  " << rmserr << std::endl;
        ti = QString ( "%1" ).arg ( rmserr );
        if ( errTab->errTabWid->item ( 1, i ) == 0 ) {
            errTab->errTabWid->setItem ( 1, i, new QTableWidgetItem ( ti ) );
        } else {
            errTab->errTabWid->item ( 1, i )->setText ( ti );
        }
        //maxval
        //::cout << "maxval =  " << maxval << std::endl;
        ti = QString ( "%1" ).arg ( maxval );
        if ( errTab->errTabWid->item ( 2, i ) == 0 ) {
            errTab->errTabWid->setItem ( 2, i, new QTableWidgetItem ( ti ) );
        } else {
            errTab->errTabWid->item ( 2, i )->setText ( ti );
        }
        //minval
        //std::cout << "minval =  " << minval << std::endl;
        ti = QString ( "%1" ).arg ( minval );
        if ( errTab->errTabWid->item ( 3, i ) == 0 ) {
            errTab->errTabWid->setItem ( 3, i, new QTableWidgetItem ( ti ) );
        } else {
            errTab->errTabWid->item ( 3, i )->setText ( ti );
        }
        //totvar
        //std::cout << "totvar =  " << totvar << std::endl;
        ti = QString ( "%1" ).arg ( totvar );
        if ( errTab->errTabWid->item ( 4, i ) == 0 ) {
            errTab->errTabWid->setItem ( 4, i, new QTableWidgetItem ( ti ) );
        } else {
            errTab->errTabWid->item ( 4, i )->setText ( ti );
        }
    }
    errTab->errTabWid->setHorizontalHeaderLabels ( colname );
}


void pde1d::setSize ( int ivalue )
{
    if ( N == ( size_t ) ivalue ) return;
    N = ivalue;
    if ( eeWidgets.empty() ) return;
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->setSize ( N );
    }

    stop = control->intTimeSteps->value();
    replot ( "Initial Condition" );
}

void pde1d::setCycles ( double value )
{
    if ( value == cycles ) return;
    cycles = value;
    double bound;
    double acyc = ( cycles > 0 ) ? cycles : -cycles;
    if ( acyc < .65 ) {
        if ( cycles < 0.5 ) {
            if ( cycles < 0.15 ) {
                if ( cycles < 0.0 ) {
                    if ( cycles < -0.15 ) {
                        if ( cycles < -0.5 ) {
                            if ( cycles < -0.65 ) {
                            } else { //  -0.65 to -0.5
                                bound = ( -0.5 - cycles ) / 0.15 + 0.2;
                                qwtPlot->setAxisScale ( 0, -1.2, bound );
                            }
                        } else { //  -0.5 to -0.15
                            qwtPlot->setAxisScale ( 0, -1.1, 0.1 );
                        }
                    } else { // -.15 to 0.0
                        bound = cycles / 0.15 - 0.1;
                        qwtPlot->setAxisScale ( 0, bound, 0.1 );
                    }
                } else { // 0 to 0.15
                    bound = cycles / 0.15 + 0.1;
                    qwtPlot->setAxisScale ( 0, -0.1, bound );
                }
            } else { // 0.15 to 0.5
                qwtPlot->setAxisScale ( 0, -0.1, 1.1 );
            }
        } else { // 0.5 to 0.65
            bound = -1.2 + ( 0.65 - cycles ) / 0.15;
            qwtPlot->setAxisScale ( 0, bound, 1.2 );
        }
    } else {
        qwtPlot->setAxisScale ( 0, -1.2, 1.2 );
    }

    if ( eeWidgets.empty() ) return;
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->initSin ( cycles );
    }
    stop = control->intTimeSteps->value();
    replot ( "Initial Condition" );
}

void pde1d::setCFL ( double value )
{
    std::cout << "pde1d::setCFL( " << value << " )\n";
    if ( value == CFL ) return;
    CFL = value;
    if ( eeWidgets.empty() ) return;
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->setCFL ( CFL );
    }
}

/*
 * void pde1d::setAlpha() {
  double alf = 1.0;
  alf = alphaLineEdit->text().toDouble();
  if(alf == lw->getAlpha() ) return;
	 lw->setAlpha(alf);
}
*/

void pde1d::setStop ( int value )
{
    if ( eeWidgets.empty() ) {
        stop = value;
    } else {
        stop = value + eeWidgets[0]->getCurrentStep();
    }
}

void pde1d::setStep ( int value )
{
    step = value;
}

void pde1d::run()
{
    if ( eeWidgets.empty() ) return;
    updateNames();
    if ( eeWidgets[0]->getCurrentStep() >= stop ) {
        metrics();
        stop = eeWidgets[0]->getCurrentStep() + control->intTimeSteps->value();
        return;
    }
    int iters = step;
    if ( eeWidgets[0]->getCurrentStep() + step > stop ) {
        iters = stop - eeWidgets[0]->getCurrentStep();
    }
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->step ( iters );
    }
    std::ostringstream s1;
    s1 << "Cycle " << std::fixed << std::setw ( 10 ) << std::setprecision ( 3 ) << ( eeWidgets[0]->getTravel() - 0.5 );
    replot ( s1.str().c_str() );
    QTimer::singleShot ( 10, this, SLOT ( run() ) );
}

void pde1d::reset()
{
    if ( eeWidgets.empty() ) return;
    updateNames();
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->initSin ( cycles );
    }
    stop = control->intTimeSteps->value();
    replot ( "Initial Condition" );
}

void pde1d::savePlot ( const QString& fileName, const char* format, int quality )
{
    std::cout << "width= " << qwtPlot->width() << std::endl;
    QImage plot ( size(), QImage::Format_RGB666 );
    render ( &plot );
    plot.save ( fileName, format, quality );
}

void pde1d::saveImage()
{
    QString fileName = QFileDialog::getSaveFileName ( this, tr ( "Save Image" ), "", tr ( "Images (*.png *.xpm *.jpg)" ) );
    if ( fileName.length() == 0 ) return;
    std::cout << fileName.toUtf8().constData() << std::endl;
    if ( fileName.contains ( QRegExp ( "//.(png|xpm|jpg)$", Qt::CaseInsensitive ) ) ) {
        savePlot ( fileName );
    } else {
        savePlot ( fileName, "PNG", 100 );
    }
}


void pde1d::addSolver ( int index )
{
  ///NOTE should change to addSolver(QString) and add combo items in constructor
    QDockWidget *qdw;
    QString qtitle;
    switch ( index ) {
    case 1:
        EEWidget * eew;
        eew = new EEWidget();
        eew->setup ( N, cycles, CFL );
        eew->setId ( dockID++ );
        connect ( eew, SIGNAL ( dockClose ( int ) ), this, SLOT ( removeSolver ( int ) ) );
        this->addDockWidget ( Qt::LeftDockWidgetArea, eew );
	this->tabifyDockWidget(control,eew);
        qdw = dynamic_cast<QDockWidget *> ( eew );
        if ( qdw == NULL ) {
            delete eew;
            return;
        }
        eeWidgets.push_back ( eew );
        break;
        //default:
    case 2:
        LeastSqrWidget *lsw;
        lsw = new LeastSqrWidget ( this );
        lsw->setup ( N, cycles, CFL );
        lsw->setId ( dockID++ );
        lsw->setAlpha ( 0.5 );
        connect ( lsw, SIGNAL ( dockClose ( int ) ), this, SLOT ( removeSolver ( int ) ) );
        this->addDockWidget ( Qt::LeftDockWidgetArea, lsw );
	this->tabifyDockWidget(control,lsw);
        qdw = dynamic_cast<QDockWidget *> ( lsw );
        if ( qdw == NULL ) {
            delete lsw;
            return;
        }
        eeWidgets.push_back ( lsw );
        break;
    case 3:
        SimpImpWidget *siw;
        siw = new SimpImpWidget ( this );
        siw->setup ( N, cycles, CFL );
        siw->setId ( dockID++ );
        siw->setImpl ( 0.5 );
        siw->setUpwind ( 0.0 );
        connect ( siw, SIGNAL ( dockClose ( int ) ), this, SLOT ( removeSolver ( int ) ) );
        this->addDockWidget ( Qt::LeftDockWidgetArea, siw );
	this->tabifyDockWidget(control,siw);
        qdw = dynamic_cast<QDockWidget *> ( siw );
        if ( qdw == NULL ) {
            delete siw;
            return;
        }
        eeWidgets.push_back ( siw );
        break;
    case 4:
        FEMWidget *femw;
        femw = new FEMWidget ( this );
        femw->setup ( N, cycles, CFL );
        femw->setId ( dockID++ );
        femw->setAlpha ( 0.5 );
        femw->setBasis ( 0 );
        connect ( femw, SIGNAL ( dockClose ( int ) ), this, SLOT ( removeSolver ( int ) ) );
        this->addDockWidget ( Qt::LeftDockWidgetArea, femw );
	this->tabifyDockWidget(control,femw);
        qdw = dynamic_cast<QDockWidget *> ( femw );
        if ( qdw == NULL ) {
            delete femw;
            return;
        }
        eeWidgets.push_back ( femw );
        break;
    case 5:
        ImpWidget *lsw2;
        lsw2 = new ImpWidget ( this );
        lsw2->setup ( N, cycles, CFL );
        lsw2->setId ( dockID++ );
        lsw2->setImplicit ( 5/12.0 );
        lsw2->setBackward ( -1/12.0 );
        lsw2->setBasis ( 0 );
        connect ( lsw2, SIGNAL ( dockClose ( int ) ), this, SLOT ( removeSolver ( int ) ) );
        this->addDockWidget ( Qt::LeftDockWidgetArea, lsw2 );
	this->tabifyDockWidget(control,lsw2);
        qdw = dynamic_cast<QDockWidget *> ( lsw2 );
        if ( qdw == NULL ) {
            delete lsw2;
            return;
        }
        eeWidgets.push_back ( lsw2 );
        break;
    case 6:
        RKWidget *rkw;
        rkw = new RKWidget ( this );
        rkw->setup ( N, cycles, CFL );
        rkw->setId ( dockID++ );
        rkw->setBasis ( 0 );
        connect ( rkw, SIGNAL ( dockClose ( int ) ), this, SLOT ( removeSolver ( int ) ) );
        this->addDockWidget ( Qt::LeftDockWidgetArea, rkw );
	this->tabifyDockWidget(control,rkw);
        qdw = dynamic_cast<QDockWidget *> ( rkw );
        if ( qdw == NULL ) {
            delete rkw;
            return;
        }
        eeWidgets.push_back ( rkw );
        break;
    default:
        return;
    }
    metrics();
    //this->addDockWidget ( Qt::LeftDockWidgetArea, qdw );
    //toolBox->addItem ( qw, qtitle );
    //connect ( qdw, SIGNAL ( destroyed ( QObject * ) ), this, SLOT ( removeSolver ( QObject * obj ) ) );
    //control->removeSolvCombo->addItem ( qtitle );
    reset();
}

void pde1d::removeSolver ( int id )
{
    size_t i;
    SolvWidget *bye;
    for ( i = 0; i < eeWidgets.size(); i++ ) {
        if ( id == eeWidgets[i]->getId() ) {
            bye = eeWidgets[i];
            eeWidgets.erase ( eeWidgets.begin() + i );
            delete bye;
            metrics();
            if ( eeWidgets.empty() ) return;
            std::ostringstream s1;
            s1 << "Cycle " << std::fixed << std::setw ( 10 ) << std::setprecision ( 3 ) << ( eeWidgets[0]->getTravel() - 0.5 );
            replot ( s1.str().c_str() );
            //QTimer::singleShot ( 10, this, SLOT ( run() ) );
            return;
        }
    }
}


void pde1d::stopRun()
{
    if ( eeWidgets.empty() ) return;
    stop = eeWidgets[0]->getCurrentStep();
}

void pde1d::updateNames()
{
    if ( eeWidgets.empty() ) return;
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        //toolBox->setItemText ( i + 1, eeWidgets[i]->getTitle() );
        //control->removeSolvCombo->setItemText ( i + 1, eeWidgets[i]->getTitle() );
    }
}
