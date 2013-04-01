#include "pde1d.h"
#include "solvwidget.h"
#include "eewidget.h"
#include "leastsqrwidget.h"
#include "simpimpwidget.h"
#include "femwidget.h"
#include "impwidget.h"
#include "rkwidget.h"
#include "envwidget.h"
#include "specwidget.h"
#include "pswidget.h"
#include "idealwidget.h"
#include "myinputs.h"

#include <QtGui/QLabel>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QAction>
#include <QtGui/QPen>
#include <QTimer>
#include <QFileDialog>
#include <QRegExp>
#include <QtGui/QImage>
#include <QtGui/QStandardItem>
#include <qwt_plot.h>
#include <qwt_legend.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_renderer.h>
#include <iostream>
#include <sstream>
#include <iomanip>

const bool pdeSolvers[4][9] =  {
    {true, true, true, true, true, true, true, true, true},
    {false, false, true, false, false, true, false, false, true},
    {false, false, true, false, false, true, false, false, true},
    {false, false, true, false, false, true, false, false, true}
};
/* "Euler Explicit"
   "Least Sqr Plus"
   "Simple Implicit"
   "Finite Element Method"
   "General Implicit"
   "Explicit Runge Kutta"
   "Envelope"
   "Spectral FFT"
   "Pseudo Spectral"
*/

/*! Main window with a qwt plot area and a control dock widget
*/
pde1d::pde1d() : QMainWindow(), Ui_MainWindow()
{
    setupUi ( this );


    dockID = 0;

    setCorner(Qt::BottomLeftCorner,Qt::LeftDockWidgetArea);

    qwtPlot->setAxisScale ( 0, -1.2, 1.2 );
    qwtPlot->setAxisScale ( 2, 0, 6.2831853 );
    qwtPlot->setCanvasBackground ( Qt::white );
    qwtPlot->setAutoDelete(false);
    
    legend = new QwtLegend();
    qwtPlot->insertLegend ( legend, QwtPlot::BottomLegend );

    errTab = new ErrTabDock ( tr ( "Metrics" ), this );
    errTab->setFeatures ( QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetMovable );
    addDockWidget ( Qt::BottomDockWidgetArea, errTab );

    control = new Controls ( tr ( "Controls" ), this );

    control->pdeBox->addItem( tr("U_t + c U_x = 0"));
    control->pdeBox->addItem( tr("U_t - d_ U_xx = 0"));
    control->pdeBox->addItem( tr("U_t+(U^2)_x=0"));
    control->pdeBox->addItem( tr("U_t+(U^2)_x-d_ U_xx"));
    connect ( control->pdeBox, SIGNAL( activated(int) ), this, SLOT( setEquation(int) ));

    connect ( control->addSolvCombo, SIGNAL ( activated ( int ) ), this, SLOT ( addSolver ( int ) ) );
    solvModel = new QStandardItemModel(this);
    control->addSolvCombo->setModel(solvModel);
    solvModel->appendRow( new QStandardItem ( tr("Add New Solver") ));
    solvModel->appendRow( new QStandardItem ( tr("Euler Explicit") ));  // addSolver(1)
    solvModel->appendRow( new QStandardItem ( tr("Least Sqr Plus") ));  // addSolver(2)
    solvModel->appendRow( new QStandardItem ( tr("Simple Implicit") )); // addSolver(3)
    solvModel->appendRow( new QStandardItem ( tr("Finite Element Method") ));// addSolver(4)
    solvModel->appendRow( new QStandardItem ( tr("General Implicit") ));// addSolver(5)
    solvModel->appendRow( new QStandardItem ( tr("Explicit Runge Kutta") ));// addSolver(6)
    solvModel->appendRow( new QStandardItem ( tr("Envelope") ));        // addSolver(7)
    solvModel->appendRow( new QStandardItem ( tr("Spectral FFT") ));    // addSolver(8)
    solvModel->appendRow( new QStandardItem ( tr("Pseudo Spectral") )); // addSolver(9)
    control->addSolvCombo->setToolTip(tr("Add a numerical solver"));


    connect ( control->viscInput, SIGNAL( valueChanged(double)), this, SLOT( setViscosity(double)));

    connect ( control->intNumberOfPoints, SIGNAL ( valueChanged ( int ) ), this, SLOT ( setSize ( int ) ) );
    control->intNumberOfPoints->setToolTip(tr("The number of points may be restricted"));
    connect ( control->cyclesInput, SIGNAL ( valueChanged ( double ) ), this, SLOT ( setCycles ( double ) ) );
    control->cyclesInput->setToolTip(tr("Positive or negative and may be fractional i.e. 0.25 - Plot scale may change"));
    connect ( control->cflInput, SIGNAL ( valueChanged ( double ) ), this, SLOT ( setCFL ( double ) ) );
    control->cflInput->setToolTip(tr("Positive = right running, many methods require |CFL| < 1.0"));
    connect ( control->intTimeSteps, SIGNAL ( valueChanged ( int ) ), this, SLOT ( setStop ( int ) ) );
    connect ( control->intPlotIncrement, SIGNAL ( valueChanged ( int ) ), this, SLOT ( setStep ( int ) ) );
    control->intPlotIncrement->setToolTip(tr("Plot incement = # points / CFL optional / cycles will stobe exact solution"));
    connect ( control->plotDelayInput, SIGNAL ( valueChanged ( int ) ), this, SLOT ( setDelay(int) ));
    control->plotDelayInput->setToolTip(tr("Extra delay (ms) to slow plot updates"));
    connect ( control->saveImageButton, SIGNAL ( clicked ( bool ) ), this, SLOT ( saveImage() ) );
    control->saveImageButton->setToolTip(tr("Saves the main window including plot to a file(default=.png)"));
    connect ( control->savePlotButton, SIGNAL ( clicked ( bool ) ), this, SLOT ( savePlot() ) );
    control->savePlotButton->setToolTip(tr("Saves the plot to a vector file(default=.svg)"));
    //connect ( control->plotDelayInput,SIGNAL( valueChanged ( int ) , this, SLOT( ..... ));
    connect ( control->resetButton, SIGNAL ( clicked ( bool ) ), this, SLOT ( reset() ) );
    control->resetButton->setToolTip(tr("Restart simulation at initial condition"));
    connect ( control->runButton, SIGNAL ( clicked ( bool ) ), this, SLOT ( run() ) );
    control->runButton->setToolTip(tr("Begin the simulation"));
    connect ( control->stopButton, SIGNAL ( clicked ( bool ) ), this, SLOT ( stopRun() ) );
    control->stopButton->setToolTip(tr("Stop the simulations at current step"));
    setTabPosition(Qt::LeftDockWidgetArea |Qt::RightDockWidgetArea, QTabWidget::West);
    control->setFeatures ( QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetMovable );
    addDockWidget ( Qt::LeftDockWidgetArea, control );
    replot ( "Initial Condition" );
    N = 100;
    cycles = 1.0;
    CFL = 1.0;
    stop = 100;
    step = 1;
    delay = -1;
    setDelay(10);
    equation = 0;
    CFL=-1;
    setCFL(1.0);
    visc = -1.0;
    setViscosity(0.001);
    
    IdealWidget * iw;
    iw = new IdealWidget(this);
    addIt(iw);
}

pde1d::~pde1d()
{}

void pde1d::replot ( const char * title )
{
    double xu[1] = {0.0};
    double uu[1] = {0.0};
    if ( eeWidgets.empty() ) return;
    std::cout << qwtPlot->autoDelete() << std::endl;
    //qwtPlot->detachItems();
    //qwtPlot->insertLegend ( legend, QwtPlot::BottomLegend );
 
    for ( size_t i = 0; i < eeWidgets.size();  i++ ) {
        if( eeWidgets[i]->canSolve(equation) ) {
            eeWidgets[i]->getCurve()->attach( qwtPlot );
        }else{
	    eeWidgets[i]->getCurve()->detach();
	}
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
    //bool notanint;
    QStringList colname;
    QTableWidgetItem oldname;
    QString ti;
    SolvWidget *solv;
    if ( eeWidgets.empty() ) return;
    if ( ( size_t ) errTab->errTabWid->columnCount() != eeWidgets.size() ) errTab->errTabWid->setColumnCount ( eeWidgets.size() );
    if ( errTab->errTabWid->rowCount() < 5 ) errTab->errTabWid->setRowCount ( 5 );
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        solv = eeWidgets[i];
        if( !(solv->canSolve(equation)) ) continue;
        //std::cout << "pde1d::metrics solver " << solv->getTitle().toLocal8Bit().data() << "  " << i << "  " << N << std::endl;

        ideal = eeWidgets[0]->getU();
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
        if( solv->isUnstable() ) {
            colname.append ( "*"+solv->getTitle()+"*" );
            ti = tr("*-*-*");
            if ( errTab->errTabWid->item ( 0, i ) == 0 ) {
                errTab->errTabWid->setItem ( 0, i, new QTableWidgetItem ( ti ) );
            } else {
                errTab->errTabWid->item ( 0, i )->setText ( ti );
            }
            if ( errTab->errTabWid->item ( 1, i ) == 0 ) {
                errTab->errTabWid->setItem ( 1, i, new QTableWidgetItem ( ti ) );
            } else {
                errTab->errTabWid->item ( 1, i )->setText ( ti );
            }
            if ( errTab->errTabWid->item ( 2, i ) == 0 ) {
                errTab->errTabWid->setItem ( 2, i, new QTableWidgetItem ( ti ) );
            } else {
                errTab->errTabWid->item ( 2, i )->setText ( ti );
            }
            if ( errTab->errTabWid->item ( 3, i ) == 0 ) {
                errTab->errTabWid->setItem ( 3, i, new QTableWidgetItem ( ti ) );
            } else {
                errTab->errTabWid->item ( 3, i )->setText ( ti );
            }
            if ( errTab->errTabWid->item ( 4, i ) == 0 ) {
                errTab->errTabWid->setItem ( 4, i, new QTableWidgetItem ( ti ) );
            } else {
                errTab->errTabWid->item ( 4, i )->setText ( ti );
            }
            continue;
        }

        colname.append ( solv->getTitle() );
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

    stop = control->intTimeSteps->getValue();
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
    stop = control->intTimeSteps->getValue();
    replot ( "Initial Condition" );
}

void pde1d::setCFL ( double value )
{
    if ( value == CFL ) return;
    CFL = value;
    if ( eeWidgets.empty() ) return;
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->setCFL ( CFL );
    }
    control->cflInput->setValue(CFL);
}

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
    int istab=-1;
    if ( eeWidgets.empty() ) return;
    for( int i=0; i< eeWidgets.size(); i++) {
        if( !(eeWidgets[i]->isUnstable()) ) {
            istab=i;
            break;
        }
    }
    if(istab == -1) return;
    //updateNames();
    if ( eeWidgets[istab]->getCurrentStep() >= stop ) {
        metrics();
        stop = eeWidgets[istab]->getCurrentStep() + control->intTimeSteps->getValue();
        return;
    }
    int iters = step;
    if ( eeWidgets[istab]->getCurrentStep() + step > stop ) {
        iters = stop - eeWidgets[istab]->getCurrentStep();
    }
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->step ( iters );
    }
    for( int i=0; i< eeWidgets.size(); i++) {
        if( !(eeWidgets[i]->isUnstable()) ) {
            istab=i;
            break;
        }
    }
    std::ostringstream s1;
    s1 << "Cycle " << std::fixed << std::setw ( 10 ) << std::setprecision ( 3 ) << ( eeWidgets[istab]->getTravel() - 0.5 );
    replot ( s1.str().c_str() );
    QTimer::singleShot ( delay, this, SLOT ( run() ) );
}

void pde1d::reset()
{
    if ( eeWidgets.empty() ) return;
    updateNames();
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->initSin ( cycles );
    }
    stop = control->intTimeSteps->getValue();
    replot ( "Initial Condition" );
    metrics();
}

void pde1d::savePlot ( )
{
    int ipt;
    QString pfile = QFileDialog::getSaveFileName ( this, tr ( "Save Plot" ), "", tr ( "Images ( *.pdf *.svg *.ps *.eps *.png )" )  );
    QString form = "svg";
    QSizeF qs(161.80,100.0);
    QwtPlotRenderer rend(this);
    if ( pfile.contains ( QRegExp ( "\\.(pdf|svg|ps)$",Qt::CaseInsensitive  ) ) ) {
        ipt = pfile.lastIndexOf(".");
        ipt = pfile.size() - ipt -1;
        form = pfile;
        form = form.right(ipt);
        rend.renderDocument(qwtPlot,pfile,form,qs);
    } else if ( pfile.contains ( QRegExp ( "\\.(bmp|jpg|jpeg|png|ppm|tiff|xbm|xpm)$", Qt::CaseInsensitive ) ) ) {
        ipt = pfile.lastIndexOf(".");
        ipt = pfile.size() - ipt -1;
        form = pfile;
        form = form.right(ipt);
        rend.renderDocument(qwtPlot,pfile,form,qs);
    } else {
        rend.renderDocument(qwtPlot,pfile,"svg",qs);
    }
}

void pde1d::saveImage()
{
    QString fileName = QFileDialog::getSaveFileName ( this, tr ( "Save Image" ), "", tr ( "Images ( *.jpg *.jpeg *.png *.ppm *.tiff *.xbm *.xpm *.bmp )" ) );
    QImage plot( size(), QImage::Format_RGB666 );
    render ( &plot );
    if ( fileName.length() == 0 ) return;
    std::cout << fileName.toUtf8().constData() << std::endl;
    if ( fileName.contains ( QRegExp ( "\\.(bmp|jpg|jpeg|png|ppm|tiff|xbm|xpm)$", Qt::CaseInsensitive ) ) ) {
        plot.save( fileName );
    } else {
        plot.save( fileName, "PNG", -1);
    }
}


void pde1d::addSolver ( int index )
{
    QString qtitle;
    switch ( index ) {
    case 1:
        EEWidget * eew;
        eew = new EEWidget();
        addIt(eew);
        break;
        //default:
    case 2:
        LeastSqrWidget *lsw;
        lsw = new LeastSqrWidget ( this );
        addIt(lsw);
        break;
    case 3:
        SimpImpWidget *siw;
        siw = new SimpImpWidget ( this );
        addIt(siw);
        break;
    case 4:
        FEMWidget *femw;
        femw = new FEMWidget ( this );
        addIt(femw);
        break;
    case 5:
        ImpWidget *lsw2;
        lsw2 = new ImpWidget ( this );
        addIt(lsw2);
        break;
    case 6:
        RKWidget *rkw;
        rkw = new RKWidget ( this );
        addIt(rkw);
        break;
    case 7:
        EnvWidget *envw;
        envw = new EnvWidget ( this );
        addIt(envw);
        break;
    case 8:
        SpecWidget *specw;
        specw= new SpecWidget ( this );
        addIt(specw);
        break;
    case 9:
        PSWidget *psw;
        psw = new PSWidget ( this );
        addIt(psw);
        break;
    default:
        return;
    }
    metrics();
    reset();
    control->addSolvCombo->setCurrentIndex(0);
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

void pde1d::setEquation(int value) {
    equation=value;
    control->pdeBox->setCurrentIndex(equation);
    for( int i=0; i<9; i++ ) {
        solvModel->item(i+1)->setEnabled(pdeSolvers[equation][i]);
    }
    if ( eeWidgets.empty() ) return;
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->setEquation ( value);
    }
}

void pde1d::setViscosity(double value) {
    if(visc == value ) return;
    visc=value;
    if ( eeWidgets.empty() ) return;
    for ( size_t i = 0; i < eeWidgets.size(); i++ ) {
        eeWidgets[i]->setViscosity ( value);
    }
    control->viscInput->setValue(visc);
}

void pde1d::addIt(SolvWidget* solver) {
    QDockWidget *qdw;
    solver->setup ( N, cycles, CFL );
    solver->setId ( dockID++ );
    solver->setEquation(equation);
    solver->setViscosity(visc);
    connect ( solver, SIGNAL ( dockClose ( int ) ), this, SLOT ( removeSolver ( int ) ) );
    this->addDockWidget ( Qt::LeftDockWidgetArea, solver );
    this->tabifyDockWidget(control,solver);
    //solver->getCurve()->attach( qwtPlot );
    qdw = dynamic_cast<QDockWidget *> ( solver );
    if ( qdw == NULL ) {
        delete solver;
        return;
    }
    eeWidgets.push_back ( solver );
    replot("Initial Conditions");
}

void pde1d::setDelay(int value) {
    if(value == delay ) return;
    delay = value;
    control->plotDelayInput->setValue(delay);
}
