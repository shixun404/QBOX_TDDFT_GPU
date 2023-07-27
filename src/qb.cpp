////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008-2020 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// qb.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

#include <sys/utsname.h>
#include <unistd.h>
#include <getopt.h>
#include <cstdlib>
#include <cassert>
#include <fstream>

#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"
#if USE_UUID
#include "uuid_str.h"
#endif
#ifdef _OPENMP
#include "omp.h"
#endif
#include "MPIdata.h"

#include "Timer.h"
#include "UserInterface.h"
#include "Sample.h"

#include "AngleCmd.h"
#include "AtomCmd.h"
#include "ComputeMLWFCmd.h"
#include "ConstraintCmd.h"
#include "DistanceCmd.h"
#include "ExtForceCmd.h"
#include "FoldInWsCmd.h"
#include "HelpCmd.h"
#include "KpointCmd.h"
#include "ListAtomsCmd.h"
#include "ListSpeciesCmd.h"
#include "LoadCmd.h"
#include "MoveCmd.h"
#include "PartialChargeCmd.h"
#include "PlotCmd.h"
#include "PrintCmd.h"
#include "QuitCmd.h"
#include "RandomizeRCmd.h"
#include "RandomizeVCmd.h"
#include "RandomizeWfCmd.h"
#include "ResetRotationCmd.h"
#include "ResetVcmCmd.h"
#include "RescaleVCmd.h"
#include "ResponseCmd.h"
#include "RseedCmd.h"
#include "RunCmd.h"
#include "SaveCmd.h"
#include "SetCmd.h"
#include "SetVelocityCmd.h"
#include "SpeciesCmd.h"
#include "SpectrumCmd.h"
#include "StatusCmd.h"
#include "StrainCmd.h"
#include "TorsionCmd.h"
#include "BisectionCmd.h"
#include "RTTDCmd.h"

#include "AlphaPBE0.h"
#include "AlphaRSH.h"
#include "AtomsDyn.h"
#include "BetaRSH.h"
#include "BlHF.h"
#include "BtHF.h"
#include "Cell.h"
#include "CellDyn.h"
#include "CellLock.h"
#include "CellMass.h"
#include "ChargeMixCoeff.h"
#include "ChargeMixNdim.h"
#include "ChargeMixRcut.h"
#include "Debug.h"
#include "Dspin.h"
#include "Ecut.h"
#include "Ecutprec.h"
#include "Ecuts.h"
#include "Efield.h"
#include "ForceTol.h"
#include "Polarization.h"
#include "Emass.h"
#include "ExtStress.h"
#include "FermiTemp.h"
#include "IterCmd.h"
#include "IterCmdPeriod.h"
#include "LockCm.h"
#include "Dt.h"
#include "MuRSH.h"
#include "Nempty.h"
#include "NetCharge.h"
#include "Nspin.h"
#include "Occ.h"
#include "RefCell.h"
#include "ScfTol.h"
#include "Stress.h"
#include "StressTol.h"
#include "Thermostat.h"
#include "ThTemp.h"
#include "ThTime.h"
#include "ThWidth.h"
#include "Vext.h"
#include "WfDiag.h"
#include "WfDyn.h"
#include "Xc.h"
#include "RTPropagator.h"
#include "RTAS.h"
#include "RTMD.h"
#include "RTVP.h"
#include "RTVpEq.h"
#include "RTVpAmp.h"
#include "RTVpFreq.h"
#include "RTVpSigma.h"
#include "RTVpGauss.h"
#include "RTVpLength.h"
#include "RTVpDelay.h"
#include "RTVpInd.h"
#include "RTDeltaKick.h"
#include "RTDeltaKickAmp.h"
#include "RTDt.h"
#include "RTScfTol.h"
#include "RTRhoTol.h"
#include "RTProj.h"
#include "RTProjType.h"
#include "RTProjBd.h"
#include "RTProjKp.h"
#include "RTAndersonCoeff.h"
#include "RTAndersonDim.h"


#if TUNER
#include "Executor.hpp"
#include "utilities.hpp"
#include "dftuningqb_util.hpp"
#include <memory>
using namespace DFTuning;
#include "FourierTransform.h"

#endif



void usage(void)
{
  cerr << " use:" << endl;
  cerr << "   qb [-nstb nstb] [-nkpb nkpb] [-nspb nspb] inputfile" << endl;
  cerr << "   qb -server [-nstb nstb] [-nkpb nkpb] [-nspb nspb] "
       << "inputfile outputfile" << endl;
}



#if TUNER

namespace DFTuningQB{

struct PerfParams {
        int nstb;
        int nkpb;
        int nspb;
        int ngb;
	int dscal_unroll;
	int dscal_tb;
	int dscal_tb_sm;
	int pair_unroll;
	int pair_tb;
	int pair_tb_sm;
        PerfParams(int st, int kp, int sp, int ranks, int dscal_u, int dscal_TB, int dscal_TB_SM, int pair_u, int pair_TB, int pair_TB_SM):nstb(st), nkpb(kp), nspb(sp),dscal_unroll(dscal_u),dscal_tb(dscal_TB),dscal_tb_sm(dscal_TB_SM),pair_unroll(pair_u),pair_tb(pair_TB),pair_tb_sm(pair_TB_SM){ngb=ranks/(st*kp*sp);}
//TODO: Return error if no exact division
};

class ExecutorQBox: public Executor
{
	private:
                static std::vector<PerfParams> params;
                static int ntasks;
                static int mype;
		static ExecutorQBox* _instance;
        public:
                ExecutorQBox(int iters):ExecutorQBox(iters,1,0){}
                ExecutorQBox(int iters,int batches):ExecutorQBox(iters, batches,0){}
                ExecutorQBox(int iters, int batches, int batch_index):Executor(iters,4,batch_index){
 			assert(batches==1);
 			std::map<int,std::vector<std::string>> p;
                        readFromCSV(setFilename(PATH_INPUT,SUFFIX_INPUT,batch_index),';',batches,p);
                        int cont=0;
                        while(cont<batches){
                                PerfParams pp(std::stoi(p[cont][1]),std::stoi(p[cont][2]),std::stoi(p[cont][3]), std::stoi(p[cont][4]),std::stoi(p[cont][5]), std::stoi(p[cont][6]),std::stoi(p[cont][7]),std::stoi(p[cont][8]),std::stoi(p[cont][9]),std::stoi(p[cont][10]));
                                params.push_back(pp);
                                cont++;
                        }
			setenv("QBOX_DSCAL_unroll",p[0][5].c_str(),1);
			setenv("QBOX_DSCAL_TB",p[0][6].c_str(),1);
			setenv("QBOX_DSCAL_TB_SM",p[0][7].c_str(),1);

			setenv("QBOX_PAIR_unroll",p[0][8].c_str(),1);
                        setenv("QBOX_PAIR_TB",p[0][9].c_str(),1);
                        setenv("QBOX_PAIR_TB_SM",p[0][10].c_str(),1);

			setenv("QBOX_ZCOPY_unroll",p[0][11].c_str(),1);
		        setenv("QBOX_ZCOPY_TB",p[0][12].c_str(),1);
                        setenv("QBOX_ZCOPY_TB_SM",p[0][13].c_str(),1);	

			
	                setenv("QBOX_VEC_unroll",p[0][14].c_str(),1);
                        setenv("QBOX_VEC_TB",p[0][15].c_str(),1);
                        setenv("QBOX_VEC_TB_SM",p[0][16].c_str(),1);

			setenv("QBOX_ZVEC_unroll",p[0][17].c_str(),1);
                        setenv("QBOX_ZVEC_TB",p[0][18].c_str(),1);
                        setenv("QBOX_ZVEC_TB_SM",p[0][19].c_str(),1);

			setenv("QBOX_NSTREAMS",p[0][20].c_str(),1);
                        setenv("QBOX_NBATCHES",p[0][21].c_str(),1);

			FourierTransform::initializeConst();
                        MPI_Init(NULL,NULL);
                        MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
                        MPI_Comm_rank(MPI_COMM_WORLD,&mype);
                                
                        //TODO: If output files exist, delete them before starting!!
                }
                ExecutorQBox(const ExecutorQBox &)=delete;
                ~ExecutorQBox(){
                        MPI_Finalize();
                }
                
                static int getNstb(const int i=0) {return params[i].nstb;}
                static int getNkpb(const int i=0) {return params[i].nkpb;}
                static int getNspb(const int i=0) {return params[i].nspb;}
                static int getNgb(const int i=0)  {return params[i].ngb;}
		
		static int getDscalU(const int i=0)  {return params[i].dscal_unroll;}
		static int getDscalTb(const int i=0)  {return params[i].dscal_tb;}
		static int getDscalTbSm(const int i=0)  {return params[i].dscal_tb_sm;}

                static int getMype(){return mype;}
		static int getNtasks (){return ntasks;}


                void writeOutput(const std::string& filename, const int index=1, const int iter=0) const{
                //TODO: MANAGE EXCEPTIONS and control returns values. What if breaks? Manage resource!!
                        if ( MPIdata::onpe0() )
                        {
					std::ofstream myfile;
                                	myfile.open (filename,std::ofstream::app);
                                	for(int i=0;i<index;i++)
						myfile << "Application.Time"<<i<<": " <<getTime(i,iter)<<"\n";
                               	 	myfile.close();
                        }       
                }

		static ExecutorQBox* Instance(int iters, int batches, int batch_index){
			if(_instance==0)
 	                      _instance = new ExecutorQBox(iters,batches,batch_index);
                return _instance;

		
		}
}; //end class
   //

int execute_test (int argc, char ** argv)
{

for(int dft_iterations=0;dft_iterations<DFTuning::Executor::getIters();dft_iterations++){	
  	Timer tm;
  	tm.start();

	auto ntasks= ExecutorQBox::getNtasks();
	auto mype = ExecutorQBox::getMype();
	int ngb = ntasks, nstb = 1, nkpb = 1, nspb = 1;

  // process command line arguments
  int ch;
  opterr = 0; // prevent getopt_long from writing on stderr
  int server_mode = 0;
  bool do_exit = false;

  // options descriptor
  static struct option longopts[] =
  {
    { "help",       no_argument,            NULL,            0  },
    { "server",     no_argument,            NULL,            1  },
    { "nstb",       required_argument,      NULL,            2  },
    { "nkpb",       required_argument,      NULL,            3  },
    { "nspb",       required_argument,      NULL,            4  },
    { NULL,         0,                      NULL,            0  }
  };
  while ( (ch = getopt_long_only(argc, argv, ":", longopts, NULL )) != -1 )
  {
    switch (ch)
    {
      case 0:
        // help
        do_exit = true;
      case 1:
        server_mode = 1;
        break;
      case 2:
        nstb = atoi(optarg);
        if ( nstb < 1 )
        {
          if ( mype == 0 )
            cerr << " nstb must be positive" << endl;
          do_exit = true;
        }
        break;
      case 3:
        nkpb = atoi(optarg);
        if ( nkpb < 1 )
        {
          if ( mype == 0 )
            cerr << " nkpb must be positive" << endl;
          do_exit = true;
        }
        break;
      case 4:
        nspb = atoi(optarg);
        if ( (nspb < 1) || (nspb > 2) )
        {
          if ( mype == 0 )
            cerr << " nspb must be 1 or 2" << endl;
          do_exit = true;
        }
        break;
      case ':':
        if ( mype == 0 )
          cerr << " missing option argument\n";
        do_exit = true;
        break;
      case '?':
        if ( mype == 0 )
          cerr << " unknown or ambiguous option\n";
        do_exit = true;
        break;
      default:
        if ( mype == 0 )
          cerr << " unknown option: " << argv[opterr] << endl;
        do_exit = true;
    }
  }

	argc -= optind;
  argv += optind;


  const int interactive = ( argc == 0 );

  nstb=ExecutorQBox::getNstb();
  if ( nstb < 1 )
        {
          if ( mype == 0 )
            cerr << " nstb must be positive" << endl;
          do_exit = true;
        }

  nkpb=ExecutorQBox::getNkpb();
  if ( nkpb < 1 )
        {
          if ( mype == 0 )
            cerr << " nkpb must be positive" << endl;
          do_exit = true;
        }

  nspb=ExecutorQBox::getNspb();
  if ( (nspb < 1) || (nspb > 2) )
        {
          if ( mype == 0 )
            cerr << " nspb must be 1 or 2" << endl;
          do_exit = true;
        }




 if ( do_exit )
  {
    if ( mype == 0 )
      usage();
    //MPI_Finalize();
    return 0;
  }




#ifdef DEBUG
  cout << " argc=" << argc << endl;
  for ( int iarg = 0; iarg < argc; ++iarg )
    cout << " argv[" << iarg << "]= " << argv[iarg] << endl;
  cout << " server_mode = " << server_mode << endl;
  cout << " interactive = " << interactive << endl;
  cout << " ntasks = " << ntasks << endl;
  cout << " nstb = " << nstb << endl;
  cout << " nkpb = " << nkpb << endl;
  cout << " nspb = " << nspb << endl;
#endif
	
 if ( ( ntasks % ( nstb * nkpb * nspb ) ) == 0 )
  {
    // nstb * nkpb * nspb divides ntasks
    ngb = ntasks / ( nstb * nkpb * nspb );
  }
  else
  {
    if ( mype == 0 )
      cerr << " nstb * nkpb * nspb does not divide ntasks evenly" << endl;
    //MPI_Finalize();
    return 0;
  }

  MPIdata::set(ngb,nstb,nkpb,nspb);

#ifdef DEBUG
  cout << MPIdata::rank() << ": ngb=" << ngb << " nstb=" << nstb
       << " nkpb=" << nkpb << " nspb=" << nspb << endl;

  cout << MPIdata::rank() << ": igb=" << MPIdata::igb()
       << " istb=" << MPIdata::istb()
       << " ikpb=" << MPIdata::ikpb()
       << " ispb=" << MPIdata::ispb() << endl;
#endif
 if ( MPIdata::onpe0() )
  {
    cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    cout << "<fpmd:simulation xmlns:fpmd=\"" << qbox_xmlns() << "\">" << endl;
#if USE_UUID
    cout << "<uuid> " << uuid_str() << " </uuid>" << endl;
#endif
    cout << "\n";
    cout << "                   ============================\n";
    cout << "                   I qbox "
         << setw(17) << left << release() << "   I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I http://qboxcode.org      I\n";
    cout << "                   ============================\n\n";
    cout << "\n";
    cout << "<release> " << release();
#ifdef TARGET
    cout << " " << TARGET;
#endif
#ifdef VERSION
    cout << " " << VERSION;
#endif
    cout << " </release>" << endl;

    // Identify executable name, checksum, size and link date
    if ( getenv("LOGNAME") != 0 )
      cout << "<user> " << getenv("LOGNAME") << " </user>" << endl;

    // Identify platform
    {
      struct utsname un;
      uname (&un);
      cout << "<sysname> " << un.sysname << " </sysname>" << endl;
      cout << "<nodename> " << un.nodename << " </nodename>" << endl;
    }

    cout << "<start_time> " << isodate() << " </start_time>" << endl;
    cout << " MPIdata::comm: " << ngb << "x" << nstb << "x"
         << nkpb << "x" << nspb << endl;
  }

  // Print list of node names
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  for ( int i = 0; i < MPI_MAX_PROCESSOR_NAME; i++ )
    processor_name[i] = '\0';
  char buf[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  PMPI_Get_processor_name(processor_name,&namelen);
// remove angle brackets from processor name for XML compatibility
  for ( int i = 0; i < MPI_MAX_PROCESSOR_NAME; i++ )
  {
    if ( processor_name[i] == '<' ) processor_name[i] = '(';
    if ( processor_name[i] == '>' ) processor_name[i] = ')';
  }

  int coords[4];
  MPI_Cart_coords(MPIdata::comm(),MPIdata::rank(),4,coords);

  if ( MPIdata::onpe0() )
  {
    cout << "<mpi_processes count=\"" << MPIdata::size() << "\">" << endl;
    cout << "<process id=\"" << MPIdata::rank() << "\"> " << processor_name
         << " </process>"
         << " (" << coords[3] << "," << coords[2]
         << "," << coords[1] << "," << coords[0] << ")" << endl;
  }
  for ( int ip = 1; ip < MPIdata::size(); ip++ )
  {
    MPI_Barrier(MPIdata::comm());
    if ( MPIdata::onpe0() )
    {
      MPI_Status status;
      MPI_Recv(&buf[0],MPI_MAX_PROCESSOR_NAME,MPI_CHAR,
                   ip,ip,MPIdata::comm(),&status);
    }
    else if ( ip == MPIdata::rank() )
    {
      // send processor name to pe0
      MPI_Send(&processor_name[0],MPI_MAX_PROCESSOR_NAME,
        MPI_CHAR,0,MPIdata::rank(),MPIdata::comm());
    }
    if ( MPIdata::onpe0() )
    {
      MPI_Status status;
      MPI_Recv(coords,4,MPI_INT,ip,ip,MPIdata::comm(),&status);
    }
    else if ( ip == MPIdata::rank() )
    {
      // send processor name to pe0
      MPI_Send(coords,4,MPI_INT,0,MPIdata::rank(),MPIdata::comm());
    }
    if ( MPIdata::onpe0() )
    {
      cout << "<process id=\"" << ip << "\"> " << buf
           << " </process>"
           << " (" << coords[3] << "," << coords[2]
           << "," << coords[1] << "," << coords[0] << ")" << endl;
    }
  }
  if ( MPIdata::onpe0() )
    cout << "</mpi_processes>" << endl;

#ifdef _OPENMP
  if ( MPIdata::onpe0() )
    cout << "<omp_max_threads> " << omp_get_max_threads()
         << " </omp_max_threads>" << endl;
#endif

  UserInterface ui;
 Sample* s = new Sample(&ui);

  ui.addCmd(new AngleCmd(s));
  ui.addCmd(new AtomCmd(s));
  ui.addCmd(new BisectionCmd(s));
  ui.addCmd(new ComputeMLWFCmd(s));
  ui.addCmd(new ConstraintCmd(s));
  ui.addCmd(new DistanceCmd(s));
  ui.addCmd(new ExtForceCmd(s));
  ui.addCmd(new FoldInWsCmd(s));
  ui.addCmd(new HelpCmd(s));
  ui.addCmd(new KpointCmd(s));
  ui.addCmd(new ListAtomsCmd(s));
  ui.addCmd(new ListSpeciesCmd(s));
  ui.addCmd(new LoadCmd(s));
  ui.addCmd(new MoveCmd(s));
  ui.addCmd(new PartialChargeCmd(s));
  ui.addCmd(new PlotCmd(s));
  ui.addCmd(new PrintCmd(s));
  ui.addCmd(new QuitCmd(s));
  ui.addCmd(new RandomizeRCmd(s));
  ui.addCmd(new RandomizeVCmd(s));
  ui.addCmd(new RandomizeWfCmd(s));
  ui.addCmd(new RescaleVCmd(s));
  ui.addCmd(new ResetRotationCmd(s));
  ui.addCmd(new ResetVcmCmd(s));
  ui.addCmd(new ResponseCmd(s));
  ui.addCmd(new RseedCmd(s));
  ui.addCmd(new RunCmd(s));
  ui.addCmd(new SaveCmd(s));
  ui.addCmd(new SetCmd(s));
  ui.addCmd(new SetVelocityCmd(s));
  ui.addCmd(new SpeciesCmd(s));
  ui.addCmd(new SpectrumCmd(s));
  ui.addCmd(new StatusCmd(s));
  ui.addCmd(new StrainCmd(s));
  ui.addCmd(new TorsionCmd(s));
  ui.addCmd(new RTTDCmd(s));

  ui.addVar(new AlphaPBE0(s));
  ui.addVar(new AlphaRSH(s));
  ui.addVar(new AtomsDyn(s));
  ui.addVar(new BetaRSH(s));
  ui.addVar(new BlHF(s));
  ui.addVar(new BtHF(s));
  ui.addVar(new Cell(s));
  ui.addVar(new CellDyn(s));
  ui.addVar(new CellLock(s));
  ui.addVar(new CellMass(s));
  ui.addVar(new ChargeMixCoeff(s));
  ui.addVar(new ChargeMixNdim(s));
  ui.addVar(new ChargeMixRcut(s));
  ui.addVar(new Debug(s));
  ui.addVar(new Dt(s));
  ui.addVar(new Ecut(s));
  ui.addVar(new Ecutprec(s));
  ui.addVar(new Ecuts(s));
  ui.addVar(new Efield(s));
  ui.addVar(new ForceTol(s));
  ui.addVar(new Polarization(s));
 ui.addVar(new Emass(s));
  ui.addVar(new ExtStress(s));
  ui.addVar(new FermiTemp(s));
  ui.addVar(new IterCmd(s));
  ui.addVar(new IterCmdPeriod(s));
  ui.addVar(new LockCm(s));
  ui.addVar(new MuRSH(s));
  ui.addVar(new Nempty(s));
  ui.addVar(new NetCharge(s));
  ui.addVar(new Nspin(s));
  ui.addVar(new Occ(s));
  ui.addVar(new Dspin(s));
  ui.addVar(new RefCell(s));
  ui.addVar(new ScfTol(s));
  ui.addVar(new Stress(s));
  ui.addVar(new StressTol(s));
  ui.addVar(new Thermostat(s));
  ui.addVar(new ThTemp(s));
  ui.addVar(new ThTime(s));
  ui.addVar(new ThWidth(s));
  ui.addVar(new Vext(s));
  ui.addVar(new WfDiag(s));
  ui.addVar(new WfDyn(s));
  ui.addVar(new Xc(s));
  ui.addVar(new RTPropagator(s));
  ui.addVar(new RTAS(s));
  ui.addVar(new RTMD(s));
  ui.addVar(new RTVP(s));
  ui.addVar(new RTVpEq(s));
  ui.addVar(new RTVpAmp(s));
  ui.addVar(new RTVpFreq(s));
  ui.addVar(new RTVpSigma(s));
  ui.addVar(new RTVpGauss(s));
  ui.addVar(new RTVpLength(s));
  ui.addVar(new RTVpDelay(s));
  ui.addVar(new RTVpInd(s));
  ui.addVar(new RTDeltaKick(s));
  ui.addVar(new RTDeltaKickAmp(s));
  ui.addVar(new RTDt(s));
  ui.addVar(new RTScfTol(s));
  ui.addVar(new RTRhoTol(s));
  ui.addVar(new RTProj(s));
  ui.addVar(new RTProjType(s));
  ui.addVar(new RTProjBd(s));
  ui.addVar(new RTProjKp(s));
  ui.addVar(new RTAndersonCoeff(s));
  ui.addVar(new RTAndersonDim(s));

if ( server_mode )
  {
    // server mode
    // input and output files expected as arguments
    if ( argc < 2 )
    {
      if ( MPIdata::onpe0() )
      {
        cout << " server mode requires two arguments" << endl;
        usage();
      }
//      MPI_Finalize();
      return 0;
    }
    string inputfilename(argv[0]);
    string outputfilename(argv[1]);
    if ( MPIdata::onpe0() )
    {
      cout << " server mode" << endl;
      cout << " input file:  " << inputfilename << endl;
      cout << " output file: " << outputfilename << endl;
    }
    bool echo = true;
    ui.processCmdsServer(inputfilename, outputfilename, "[qbox]", echo);
  }
  else if ( interactive )
  {
    // interactive mode
    assert(argc==0);
    // use standard input
    bool echo = !isatty(0);
    ui.processCmds(cin, "[qbox]", echo);
  }
  else
  {
    // cmd line: qb inputfilename
    // input file expected as a command line argument
    assert(argc >= 1);
    bool echo = true;
    ifstream in;
    int file_ok = 0;
    if ( MPIdata::onpe0() )
    {
      in.open(argv[0],ios::in);
      if ( in )
      {
        // file was opened on process 0
        file_ok = 1;
      }
    }
    MPI_Bcast(&file_ok,1,MPI_INT,0,MPI_COMM_WORLD);
if ( file_ok )
    {
      ui.processCmds(in, "[qbox]", echo);
    }
    else
    {
      if ( MPIdata::onpe0() )
        cout << " Could not open input file " << argv[0] << endl;
  //    MPI_Finalize();
      return 0;
    }
  }

  // exit using the quit command when processCmds returns
  Cmd *c = ui.findCmd("quit");
  c->action(1,NULL);

  if ( MPIdata::onpe0() )
  {
    cout << "<real_time> " << tm.real() << " </real_time>" << endl;
    cout << "<end_time> " << isodate() << " </end_time>" << endl;
    cout << "</fpmd:simulation>" << endl;
    DFTuning::Executor::setTime(tm.real(),3,dft_iterations);
  }

  delete s;

//  MPI_Finalize();


}//end iters

 return 0;

}//end execution method	




} //end namespace


const std::shared_ptr<Executor> factoryFun(const int iters, const int batches, const int batch_index){

        std::shared_ptr<Executor> p(DFTuningQB::ExecutorQBox::Instance(iters,batches,batch_index));
        return p;
}

DFTuningQB::ExecutorQBox* DFTuningQB::ExecutorQBox::_instance=0;
std::vector<DFTuningQB::PerfParams> DFTuningQB::ExecutorQBox::params;
std::vector<std::vector<double>> DFTuning::Executor::time;
int DFTuning::Executor::iters_=-1;
int DFTuning::Executor::batch_index_=-1;
int DFTuning::Executor::batches_=-1;
int DFTuningQB::ExecutorQBox::ntasks=-1;
int DFTuningQB::ExecutorQBox::mype=-1;


int FourierTransform::nstreams=-1;
int FourierTransform::nbatches=-1;

int main(int argc, char **argv, char **envp)
{
        int niters= atoi(getEnv("DFTUNING_NITERS"));
        int argcc=0;
        argcc=argc-1;
        int batch_index=-1;
        batch_index=atoi(argv[argcc]);

        auto o = factoryFun(niters,1,batch_index);
	o->execution(DFTuningQB::execute_test,argcc,argv);
         for(int i=0;i<(o->getIters());i++){
                 o->writeOutput(setFilename(DFTuningQB::PATH_OUTPUT,DFTuningQB::SUFFIX_OUTPUT,batch_index),4,i);
         }
         return 0;
}




#else



int main(int argc, char **argv, char **envp)
{
  Timer tm;
  tm.start();

  MPI_Init(&argc,&argv);

  int ntasks;
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  int mype;
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  // default values for number of blocks
  // ngb: number of G vector blocks
  // nstb: number of states blocks
  // nkpb: number of kpoint blocks
  // nspb: number of spin blocks
  int ngb = ntasks, nstb = 1, nkpb = 1, nspb = 1;

  // process command line arguments
  int ch;
  opterr = 0; // prevent getopt_long from writing on stderr
  int server_mode = 0;
  bool do_exit = false;

  // options descriptor
  static struct option longopts[] =
  {
    { "help",       no_argument,            NULL,            0  },
    { "server",     no_argument,            NULL,            1  },
    { "nstb",       required_argument,      NULL,            2  },
    { "nkpb",       required_argument,      NULL,            3  },
    { "nspb",       required_argument,      NULL,            4  },
    { NULL,         0,                      NULL,            0  }
  };

  while ( (ch = getopt_long_only(argc, argv, ":", longopts, NULL )) != -1 )
  {
    switch (ch)
    {
      case 0:
        // help
        do_exit = true;
      case 1:
        server_mode = 1;
        break;
      case 2:
        nstb = atoi(optarg);
        if ( nstb < 1 )
        {
          if ( mype == 0 )
            cerr << " nstb must be positive" << endl;
          do_exit = true;
        }
        break;
      case 3:
        nkpb = atoi(optarg);
        if ( nkpb < 1 )
        {
          if ( mype == 0 )
            cerr << " nkpb must be positive" << endl;
          do_exit = true;
        }
        break;
      case 4:
        nspb = atoi(optarg);
        if ( (nspb < 1) || (nspb > 2) )
        {
          if ( mype == 0 )
            cerr << " nspb must be 1 or 2" << endl;
          do_exit = true;
        }
        break;
      case ':':
        if ( mype == 0 )
          cerr << " missing option argument\n";
        do_exit = true;
        break;
      case '?':
        if ( mype == 0 )
          cerr << " unknown or ambiguous option\n";
        do_exit = true;
        break;
      default:
        if ( mype == 0 )
          cerr << " unknown option: " << argv[opterr] << endl;
        do_exit = true;
    }
  }
  argc -= optind;
  argv += optind;

  if ( do_exit )
  {
    if ( mype == 0 )
      usage();
    MPI_Finalize();
    return 0;
  }

  const int interactive = ( argc == 0 );

#ifdef DEBUG
  cout << " argc=" << argc << endl;
  for ( int iarg = 0; iarg < argc; ++iarg )
    cout << " argv[" << iarg << "]= " << argv[iarg] << endl;
  cout << " server_mode = " << server_mode << endl;
  cout << " interactive = " << interactive << endl;
  cout << " ntasks = " << ntasks << endl;
  cout << " nstb = " << nstb << endl;
  cout << " nkpb = " << nkpb << endl;
  cout << " nspb = " << nspb << endl;
#endif

  // adjust ngb to satisfy ngb*nstb*nkpb*nspb == ntasks

  if ( ( ntasks % ( nstb * nkpb * nspb ) ) == 0 )
  {
    // nstb * nkpb * nspb divides ntasks
    ngb = ntasks / ( nstb * nkpb * nspb );
  }
  else
  {
    if ( mype == 0 )
      cerr << " nstb * nkpb * nspb does not divide ntasks evenly" << endl;
    MPI_Finalize();
    return 0;
  }

  MPIdata::set(ngb,nstb,nkpb,nspb);

#ifdef DEBUG
  cout << MPIdata::rank() << ": ngb=" << ngb << " nstb=" << nstb
       << " nkpb=" << nkpb << " nspb=" << nspb << endl;

  cout << MPIdata::rank() << ": igb=" << MPIdata::igb()
       << " istb=" << MPIdata::istb()
       << " ikpb=" << MPIdata::ikpb()
       << " ispb=" << MPIdata::ispb() << endl;
#endif

  if ( MPIdata::onpe0() )
  {
    cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    cout << "<fpmd:simulation xmlns:fpmd=\"" << qbox_xmlns() << "\">" << endl;
#if USE_UUID
    cout << "<uuid> " << uuid_str() << " </uuid>" << endl;
#endif
    cout << "\n";
    cout << "                   ============================\n";
    cout << "                   I qbox "
         << setw(17) << left << release() << "   I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I                          I\n";
    cout << "                   I http://qboxcode.org      I\n";
    cout << "                   ============================\n\n";
    cout << "\n";
    cout << "<release> " << release();
#ifdef TARGET
    cout << " " << TARGET;
#endif
#ifdef VERSION
    cout << " " << VERSION;
#endif
    cout << " </release>" << endl;

    // Identify executable name, checksum, size and link date
    if ( getenv("LOGNAME") != 0 )
      cout << "<user> " << getenv("LOGNAME") << " </user>" << endl;

    // Identify platform
    {
      struct utsname un;
      uname (&un);
      cout << "<sysname> " << un.sysname << " </sysname>" << endl;
      cout << "<nodename> " << un.nodename << " </nodename>" << endl;
    }

    cout << "<start_time> " << isodate() << " </start_time>" << endl;
    cout << " MPIdata::comm: " << ngb << "x" << nstb << "x"
         << nkpb << "x" << nspb << endl;
  }

  // Print list of node names
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  for ( int i = 0; i < MPI_MAX_PROCESSOR_NAME; i++ )
    processor_name[i] = '\0';
  char buf[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  PMPI_Get_processor_name(processor_name,&namelen);
  // remove angle brackets from processor name for XML compatibility
  for ( int i = 0; i < MPI_MAX_PROCESSOR_NAME; i++ )
  {
    if ( processor_name[i] == '<' ) processor_name[i] = '(';
    if ( processor_name[i] == '>' ) processor_name[i] = ')';
  }

  int coords[4];
  MPI_Cart_coords(MPIdata::comm(),MPIdata::rank(),4,coords);

  if ( MPIdata::onpe0() )
  {
    cout << "<mpi_processes count=\"" << MPIdata::size() << "\">" << endl;
    cout << "<process id=\"" << MPIdata::rank() << "\"> " << processor_name
         << " </process>"
         << " (" << coords[3] << "," << coords[2]
         << "," << coords[1] << "," << coords[0] << ")" << endl;
  }
  for ( int ip = 1; ip < MPIdata::size(); ip++ )
  {
    MPI_Barrier(MPIdata::comm());
    if ( MPIdata::onpe0() )
    {
      MPI_Status status;
      MPI_Recv(&buf[0],MPI_MAX_PROCESSOR_NAME,MPI_CHAR,
                   ip,ip,MPIdata::comm(),&status);
    }
    else if ( ip == MPIdata::rank() )
    {
      // send processor name to pe0
      MPI_Send(&processor_name[0],MPI_MAX_PROCESSOR_NAME,
        MPI_CHAR,0,MPIdata::rank(),MPIdata::comm());
    }
    if ( MPIdata::onpe0() )
    {
      MPI_Status status;
      MPI_Recv(coords,4,MPI_INT,ip,ip,MPIdata::comm(),&status);
    }
    else if ( ip == MPIdata::rank() )
    {
      // send processor name to pe0
      MPI_Send(coords,4,MPI_INT,0,MPIdata::rank(),MPIdata::comm());
    }
    if ( MPIdata::onpe0() )
    {
      cout << "<process id=\"" << ip << "\"> " << buf
           << " </process>"
           << " (" << coords[3] << "," << coords[2]
           << "," << coords[1] << "," << coords[0] << ")" << endl;
    }
  }
  if ( MPIdata::onpe0() )
    cout << "</mpi_processes>" << endl;

#ifdef _OPENMP
  if ( MPIdata::onpe0() )
    cout << "<omp_max_threads> " << omp_get_max_threads()
         << " </omp_max_threads>" << endl;
#endif

  UserInterface ui;
  Sample* s = new Sample(&ui);

  ui.addCmd(new AngleCmd(s));
  ui.addCmd(new AtomCmd(s));
  ui.addCmd(new BisectionCmd(s));
  ui.addCmd(new ComputeMLWFCmd(s));
  ui.addCmd(new ConstraintCmd(s));
  ui.addCmd(new DistanceCmd(s));
  ui.addCmd(new ExtForceCmd(s));
  ui.addCmd(new FoldInWsCmd(s));
  ui.addCmd(new HelpCmd(s));
  ui.addCmd(new KpointCmd(s));
  ui.addCmd(new ListAtomsCmd(s));
  ui.addCmd(new ListSpeciesCmd(s));
  ui.addCmd(new LoadCmd(s));
  ui.addCmd(new MoveCmd(s));
  ui.addCmd(new PartialChargeCmd(s));
  ui.addCmd(new PlotCmd(s));
  ui.addCmd(new PrintCmd(s));
  ui.addCmd(new QuitCmd(s));
  ui.addCmd(new RandomizeRCmd(s));
  ui.addCmd(new RandomizeVCmd(s));
  ui.addCmd(new RandomizeWfCmd(s));
  ui.addCmd(new RescaleVCmd(s));
  ui.addCmd(new ResetRotationCmd(s));
  ui.addCmd(new ResetVcmCmd(s));
  ui.addCmd(new ResponseCmd(s));
  ui.addCmd(new RseedCmd(s));
  ui.addCmd(new RunCmd(s));
  ui.addCmd(new SaveCmd(s));
  ui.addCmd(new SetCmd(s));
  ui.addCmd(new SetVelocityCmd(s));
  ui.addCmd(new SpeciesCmd(s));
  ui.addCmd(new SpectrumCmd(s));
  ui.addCmd(new StatusCmd(s));
  ui.addCmd(new StrainCmd(s));
  ui.addCmd(new TorsionCmd(s));
  ui.addCmd(new RTTDCmd(s));

  ui.addVar(new AlphaPBE0(s));
  ui.addVar(new AlphaRSH(s));
  ui.addVar(new AtomsDyn(s));
  ui.addVar(new BetaRSH(s));
  ui.addVar(new BlHF(s));
  ui.addVar(new BtHF(s));
  ui.addVar(new Cell(s));
  ui.addVar(new CellDyn(s));
  ui.addVar(new CellLock(s));
  ui.addVar(new CellMass(s));
  ui.addVar(new ChargeMixCoeff(s));
  ui.addVar(new ChargeMixNdim(s));
  ui.addVar(new ChargeMixRcut(s));
  ui.addVar(new Debug(s));
  ui.addVar(new Dt(s));
  ui.addVar(new Ecut(s));
  ui.addVar(new Ecutprec(s));
  ui.addVar(new Ecuts(s));
  ui.addVar(new Efield(s));
  ui.addVar(new ForceTol(s));
  ui.addVar(new Polarization(s));
  ui.addVar(new Emass(s));
  ui.addVar(new ExtStress(s));
  ui.addVar(new FermiTemp(s));
  ui.addVar(new IterCmd(s));
  ui.addVar(new IterCmdPeriod(s));
  ui.addVar(new LockCm(s));
  ui.addVar(new MuRSH(s));
  ui.addVar(new Nempty(s));
  ui.addVar(new NetCharge(s));
  ui.addVar(new Nspin(s));
  ui.addVar(new Occ(s));
  ui.addVar(new Dspin(s));
  ui.addVar(new RefCell(s));
  ui.addVar(new ScfTol(s));
  ui.addVar(new Stress(s));
  ui.addVar(new StressTol(s));
  ui.addVar(new Thermostat(s));
  ui.addVar(new ThTemp(s));
  ui.addVar(new ThTime(s));
  ui.addVar(new ThWidth(s));
  ui.addVar(new Vext(s));
  ui.addVar(new WfDiag(s));
  ui.addVar(new WfDyn(s));
  ui.addVar(new Xc(s));
  ui.addVar(new RTPropagator(s));
  ui.addVar(new RTAS(s));
  ui.addVar(new RTMD(s));
  ui.addVar(new RTVP(s));
  ui.addVar(new RTVpEq(s));
  ui.addVar(new RTVpAmp(s));
  ui.addVar(new RTVpFreq(s));
  ui.addVar(new RTVpSigma(s));
  ui.addVar(new RTVpGauss(s));
  ui.addVar(new RTVpLength(s));
  ui.addVar(new RTVpDelay(s));
  ui.addVar(new RTVpInd(s));
  ui.addVar(new RTDeltaKick(s));
  ui.addVar(new RTDeltaKickAmp(s));
  ui.addVar(new RTDt(s));
  ui.addVar(new RTScfTol(s));
  ui.addVar(new RTRhoTol(s));
  ui.addVar(new RTProj(s));
  ui.addVar(new RTProjType(s));
  ui.addVar(new RTProjBd(s));
  ui.addVar(new RTProjKp(s));
  ui.addVar(new RTAndersonCoeff(s));
  ui.addVar(new RTAndersonDim(s));

  if ( server_mode )
  {
    // server mode
    // input and output files expected as arguments
    if ( argc < 2 )
    {
      if ( MPIdata::onpe0() )
      {
        cout << " server mode requires two arguments" << endl;
        usage();
      }
      MPI_Finalize();
      return 0;
    }
    string inputfilename(argv[0]);
    string outputfilename(argv[1]);
    if ( MPIdata::onpe0() )
    {
      cout << " server mode" << endl;
      cout << " input file:  " << inputfilename << endl;
      cout << " output file: " << outputfilename << endl;
    }
    bool echo = true;
    ui.processCmdsServer(inputfilename, outputfilename, "[qbox]", echo);
  }
  else if ( interactive )
  {
    // interactive mode
    assert(argc==0);
    // use standard input
    bool echo = !isatty(0);
    ui.processCmds(cin, "[qbox]", echo);
  }
  else
  {
    // cmd line: qb inputfilename
    // input file expected as a command line argument
    assert(argc >= 1);
    bool echo = true;
    ifstream in;
    int file_ok = 0;
    if ( MPIdata::onpe0() )
    {
      in.open(argv[0],ios::in);
      if ( in )
      {
        // file was opened on process 0
        file_ok = 1;
      }
    }
    MPI_Bcast(&file_ok,1,MPI_INT,0,MPI_COMM_WORLD);

    if ( file_ok )
    {
      ui.processCmds(in, "[qbox]", echo);
    }
    else
    {
      if ( MPIdata::onpe0() )
        cout << " Could not open input file " << argv[0] << endl;
      MPI_Finalize();
      return 0;
    }
  }

  // exit using the quit command when processCmds returns
  Cmd *c = ui.findCmd("quit");
  c->action(1,NULL);

  if ( MPIdata::onpe0() )
  {
    cout << "<real_time> " << tm.real() << " </real_time>" << endl;
    cout << "<end_time> " << isodate() << " </end_time>" << endl;
    cout << "</fpmd:simulation>" << endl;
  }

  delete s;

  MPI_Finalize();

  return 0;
}
#endif
