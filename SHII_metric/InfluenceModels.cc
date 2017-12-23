#include "InfluenceModels.h"

InfluenceModels::InfluenceModels() {
}

InfluenceModels::~InfluenceModels() {
	string pid = intToStr(unsigned(getpid()));
	string outfile = "temp/tmp_" + pid + ".txt";
	
	string command = string("rm -f ") + outfile ;
	system(command.c_str());
	
	//delete input;
}

void InfluenceModels::doAll(int argc, char* argv[]) 
{ 
	// Lets harcode oneAction option for now
	actionsProcessTillNow = 0;

	time (&start);
	wastedTime = 0;

	AnyOption* opt = readOptions(argc, argv);


	cout << "In function InfluenceModels::doAll " << endl;

	cout << "PID of the program: " << getpid() << endl;

	string command = string("mkdir -p ") + outdir ;
	system(command.c_str());

	command = string("mkdir -p temp") ;
	system(command.c_str());
	
//	mkdir(outdir, 0777);

    MC* mc = new MC(opt);
    mc->doAll();
    delete mc;
	

	delete opt;

	cout << "Completed doAll in InfluenceModels" << endl;
//	writeInfMatrix();

	cout << "All memory released, current usage : " << getCurrentMemoryUsage() << endl;
}


AnyOption* InfluenceModels::readOptions(int argc, char* argv[] ) 
{
	// read the command line options
	AnyOption *opt = new AnyOption();

	// ignore POSIX style options
	opt->noPOSIX(); 
	opt->setVerbose(); /* print warnings about unknown options */
	opt->autoUsagePrint(true); /* print usage for bad options */

	/* 3. SET THE USAGE/HELP   */
	opt->addUsage( "" );
	opt->addUsage( "Usage: " );
	opt->addUsage( "" );
	opt->addUsage( " -help                            Prints this help " );
	opt->addUsage( " -c <config_file>                 Specify config file " );
	opt->addUsage( " -outdir <output_dir>             Output directory ");
	opt->addUsage( " -propModel <propagation model>   IC or LT ");
	opt->addUsage( " -mcruns <MC runs>                Number of Monte Carlo runs ");
	opt->addUsage( " -celfPlus <0 or 1>               0 for CELF and 1 for CELF++ ");
	opt->addUsage( " -probGraphFile <0 or 1>          Input file for influence network ");
//	opt->addUsage( " -seedFileName <Seed Set File Name if th>   Number of Monte Carlo runs ");
	opt->addUsage( "" );

	/* 4. SET THE OPTION STRINGS/CHARACTERS */
	/* by default all  options  will be checked on the command line and from option/resource file */
	opt->setOption("propModel");
	opt->setOption("outdir");
    opt->setOption("datasetFile");

    opt->setOption("randSeedIter");
//	opt->setOption("seedFileName");
	opt->setOption("mcruns");
    opt->setOption("celfPlus");
    opt->setOption("startIt");

	/* for options that will be checked only on the command line not in option/resource file */
	opt->setCommandFlag("help");
	opt->setCommandOption("c");

	/* for options that will be checked only from the option/resource file */
	opt->setOption("probGraphFile");
	opt->setOption("labelFile");
    opt->setOption("seedFile");
	opt->setOption("seedRatio");
	opt->setOption("edgeWeight");
	// for static or dynamic friendship graph setting

	// option for dataset size
	// 0 : movielens (very small)
	// 1 : yahoo movies (intermediate)
	// 2 : extremely large (flickr)


	/* go through the command line and get the options  */
	opt->processCommandArgs( argc, argv );

	/* 6. GET THE VALUES */
	if(opt->getFlag( "help" )) {
		opt->printUsage();
		delete opt;
		exit(0);
	}

	const char* configFile = opt->getValue("c");
	if (configFile == NULL) {
		cout << "Config file not mentioned" << endl;
		opt->printUsage();
		delete opt;
		exit(0);
	}

	cout << "Config file : " << configFile << endl;

	opt->processFile(configFile);
	opt->processCommandArgs( argc, argv );


	//training_dir = opt->getValue("training_dir");
	outdir = opt->getValue("outdir");

	//	cout << "doGenuineLeaders: " << doGenuineLeaders << endl;
	//	cout << "doTribeLeaders: " << doTribeLeaders << endl;
	//cout << "training directory: " << training_dir << endl;
	cout << "output directory: " << outdir << endl;
	cout << endl << endl;
	return opt;

}
	


int main(int argc, char* argv[]) {	
	InfluenceModels *L = new InfluenceModels();
	L->doAll(argc, argv);
	//L->writeSummary(0);
	cout << "in main again" << endl;
	delete L;
	return 0;
}


