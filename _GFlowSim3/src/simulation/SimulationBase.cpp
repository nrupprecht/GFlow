#include "SimulationBase.hpp"
#include "../integrators/VelocityVerletIntegrator.hpp"
#include "../forces/ViscousDrag.hpp"
#include "../forces/TemperatureForce.hpp"

namespace GFlow {

  SimulationBase::SimulationBase() : integrator(nullptr), simData(nullptr), dataRecord(nullptr), nowrite(false), print(false), quiet(false), saveFile(""), tag(0) {};

  SimulationBase::~SimulationBase() {
    if (simData)    delete simData;
    if (dataRecord) delete dataRecord;
    if (integrator) delete integrator;
  }

  void SimulationBase::setUp(int ac, char** av) {
    argc = ac;
    argv = av;
    // Set up arg parser, file parser
    fileParser.set(argc, argv);
    parser.set(argc, argv);
    fileParser.setArgParse(&parser);
    fileParser.setDataRecord(dataRecord);
  }

  void SimulationBase::updateFlagCheck(string flag) {
    // This will set the checked flag to true
    parser.find(flag);
  }

  void SimulationBase::standardParsing(bool check) {
    // Files
    string config = "";
    string writeDirectory = "";
    string loadFile = "";
    // Record options
    bool animate = false;
    RealType lowerSizeLimit = 0;
    int recOption = 1;
    bool recPerf = false;
    bool recMvRatio = false;
    bool recDt = false;
    bool recDelay = false;
    bool recDistance = false;
    bool recBulk = false;
    bool recBulkOutline = false;
    bool recDisplacementField = false;
    bool recDisplacementColumns = false;
    bool displacementSnapshot = false;
    bool recPressField = false;
    bool recVortexField = false;
    bool trackDisplacement = false;
    RealType fps = 15;
    RealType dataFrequency = 15;
    // Stat function options
    bool maxRatio = false;
    bool ratio = false;
    bool KE = false;
    bool maxV = false;
    bool aveV = false;
    bool aveF = false;
    bool maxF = false;
    bool trackX = false;
    bool trackY = false;
    bool mixing = false;
    // Stat plot options
    bool plotVelocity = false;
    bool plotCorrelation = false;
    bool plotDensity = false;
    bool plotPressureVDepth = false;
    bool plotAlignment = false;
    // Centering
    bool center = false;
    bool cv = false;
    bool insert = false;
    // Performance options
    int adjust = -1;       // Whether to auto-adjust the time step
    int adjustDelay = -1;  // Whether to auto-adjust the update delay
    RealType maxDt = -1;
    RealType minDt = -1;
    RealType dt = -1;
    // Temperature
    RealType temperature = -1.; // The temperature of the simulation (if positive)
    // Seed rand
    unsigned int seed = 0;

    /** Check standard options **/
    // Check if we need to load a file
    parser.get("config", config);
    parser.get("writeDirectory", writeDirectory);
    parser.get("loadFile", loadFile);
    parser.get("saveFile", saveFile);
    // Animation options
    parser.get("animate", animate);
    parser.get("lowerSizeLimit", lowerSizeLimit);
    parser.get("recOption", recOption);
    parser.get("recPerf", recPerf);
    parser.get("recMvRatio", recMvRatio);
    parser.get("recDt", recDt);
    parser.get("recDelay", recDelay);
    parser.get("recDistance", recDistance);
    parser.get("recBulk", recBulk);
    parser.get("recBulkOutline", recBulkOutline);
    parser.get("recDisplacementField", recDisplacementField);
    parser.get("recDisplacementColumns", recDisplacementColumns);
    parser.get("displacementSnapshot", displacementSnapshot);
    parser.get("recPressField", recPressField);
    parser.get("recVortexField", recVortexField);
    parser.get("trackDisplacement", trackDisplacement);
    parser.get("fps", fps);
    parser.get("dataFrequency", dataFrequency);
    parser.get("maxRatio", maxRatio);
    parser.get("ratio", ratio);
    parser.get("KE", KE);
    parser.get("maxV", maxV);
    parser.get("aveV", aveV);
    parser.get("maxF", maxF);
    parser.get("aveF", aveF);
    parser.get("trackX", trackX);
    parser.get("mixing", mixing);
    parser.get("trackY", trackY);
    parser.get("plotVelocity", plotVelocity);
    parser.get("plotCorrelation", plotCorrelation);
    parser.get("plotDensity", plotDensity);
    parser.get("plotPressureVDepth", plotPressureVDepth);
    parser.get("plotAlignment", plotAlignment);
    parser.get("center", center);
    parser.get("print", print);
    parser.get("quiet", quiet);
    parser.get("nowrite", nowrite);
    // Performance options
    parser.get("adjust", adjust);
    parser.get("adjustDelay", adjustDelay);
    parser.get("maxDt", maxDt);
    parser.get("minDt", minDt);
    parser.get("dt", dt);
    // Temperature
    parser.get("temperature", temperature);
    // Tag
    parser.get("tag", tag);
    // Randomness
    parser.get("seed", seed);

    // Seed random number generator
    if (seed==0) seed = std::chrono::system_clock::now().time_since_epoch().count();
    fileParser.setSeed(seed);
    srand48( seed );
    seedNormalDistribution();

    // Set up configurations
    if (config!="") {
      try {
	fileParser.parse(config, simData, integrator);
      }
      catch (FileParser::FileDoesNotExist file) {
	if (!quiet) std::cerr << "File [" << file.name << "] does not exist. Trying [samples/" << file.name << "].\n";
	try {
	  fileParser.parse("samples/"+config, simData, integrator);
	}
	catch (FileParser::FileDoesNotExist file) {
	  std::cerr << "File [" << file.name << "] also does not exist. Exiting.\n";
	  exit(0);
	}
	if (!quiet) cout << "You're lucky, [samples/" << file.name << "] does exist. Assuming you meant that file.\n";
      }
      catch (FileParser::UnrecognizedToken token) {
	std::cerr << "Unrecognized token [" << token.token << "]. Exiting.\n";
	exit(0);
      }
      if (simData==nullptr) {
	std::cerr << "Error occured while parsing, sim data is null. Exiting.\n";
	exit(0);
      }
    }
    // Load a saved state
    else if (loadFile!="") {
      try {
	fileParser.loadFromFile(loadFile, simData, integrator);
      }
      catch (FileParser::FileDoesNotExist file) {
	std::cerr << "File [" << file.name << "] does not exist. Trying [samples/" << file.name << "]. Exiting.\n";
	exit(0);
      }
      catch (GFlow::BadIstreamRead read) {
	std::cerr << "Parsing failed when it attempted to read in a [" << read.type << "], expected [" << printChar(read.expected) << "], found [" << printChar(read.unexpected) << "]. Exiting.\n";
	exit(0);
      }
    }
    // Otherwise, load the default configuration
    else {
      try {
	fileParser.parse("samples/box.cfg", simData, integrator);
      }
      catch (FileParser::FileDoesNotExist file) {
	std::cerr << "Default file [samples/box.cfg] does not exist. Try loading a specifi\
c configuration file. Exiting.\n";
	exit(0);
      }
      catch (FileParser::UnrecognizedToken token) {
	std::cerr << "Unrecognized token [" << token.token << "]. Exiting.\n";
	exit(0);
      }
    }

    // Make sure simData is non-null
    if (simData==nullptr) {
      std::cerr << "SimData is null. Exiting." << endl;
      exit(0);
    }
    
    // Set integrator options
    if (maxDt>0) {
      VelocityVerletIntegrator* integ = nullptr;
      integ = reinterpret_cast<VelocityVerletIntegrator*>(integrator);
      if (integrator) integ->setMaxTimeStep(maxDt);
    }
    if (minDt>0) {
      VelocityVerletIntegrator* integ = nullptr;
      integ = reinterpret_cast<VelocityVerletIntegrator*>(integrator);
      if (integ) integ->setMinTimeStep(minDt);
    }
    if (dt>0) integrator->setDt(dt);

    // Temperature
    if (temperature>0) {
      simData->addExternalForce(new TemperatureForce(temperature));
      simData->addExternalForce(new ViscousDrag);
    }

    // Set up the data recorder
    if (writeDirectory!="") dataRecord->setWriteDirectory(writeDirectory);
    integrator->setDataRecord(dataRecord);
    // Set record options
    if (dataRecord) {
      // Set stat function options
      if (KE)       dataRecord->addStatFunction(StatFunc_AveKE, "KE");
      if (maxRatio) dataRecord->addStatFunction(StatFunc_MaxVelocitySigmaRatio, "maxRatio");
      if (ratio)    dataRecord->addStatFunction(StatFunc_MinSigmaVelocityRatio, "ratio");
      if (maxV)     dataRecord->addStatFunction(StatFunc_MaxSpeed, "maxV");
      if (aveV)     dataRecord->addStatFunction(StatFunc_AveSpeed, "aveV");
      if (maxF)     dataRecord->addStatFunction(StatFunc_MaxForce, "maxF");
      if (aveF)     dataRecord->addStatFunction(StatFunc_AveForce, "aveF");
      if (trackX)   dataRecord->addStatFunction(StatFunc_MaxR_PosX, "trackX");
      if (trackY)   dataRecord->addStatFunction(StatFunc_MaxR_PosY, "trackY");
      if (mixing)   dataRecord->addStatFunction(StatFunc_MixingParameter, "mixing");
      // Set stat plot options
      if (plotVelocity) dataRecord->addStatPlot(StatPlot_Velocity, RPair(0,3), 100, "velocityPlot");
      if (plotCorrelation) dataRecord->addStatPlot(StatPlot_RadialCorrelation, RPair(0,1), 100, "correlation");
      if (plotDensity) dataRecord->addStatPlot(StatPlot_DensityVsDepth, RPair(0,0), 100, "densityVsDepth");
      if (plotPressureVDepth) dataRecord->addStatPlot(StatPlot_PressureVsDepth, RPair(0,0), 100, "pressureVsDepth");
      if (plotAlignment) dataRecord->addStatPlot(StatPlot_Alignment, RPair(0,1), 100, "alignment");
      // Set recording options
      dataRecord->setRecPos(animate);
      dataRecord->setLowerSizeLimit(lowerSizeLimit);
      dataRecord->setRecOption(recOption);
      dataRecord->setRecPerf(recPerf);
      dataRecord->setRecMvRatio(recMvRatio);
      dataRecord->setRecDt(recDt);
      dataRecord->setRecDelay(recDelay);
      dataRecord->setRecDistance(recDistance);
      dataRecord->setRecBulk(recBulk);
      dataRecord->setRecBulkOutline(recBulkOutline);
      dataRecord->setRecDisplacementField(recDisplacementField);
      dataRecord->setRecDisplacementColumns(recDisplacementColumns);
      dataRecord->setDisplacementSnapshot(displacementSnapshot);
      dataRecord->setRecPressField(recPressField);
      dataRecord->setRecVortex(recVortexField);
      dataRecord->setTrackDisplacement(trackDisplacement);
      dataRecord->setAnimationDelay(1./fps);
      dataRecord->setDelay(1./dataFrequency);
      dataRecord->setCenter(center);
    }

    // Set SimData initial positions
    simData->setInitialPositions();

    // Set run time
    if (adjust>-1) {
      auto v = dynamic_cast<VelocityVerletIntegrator*>(integrator);
      if (v) v->setAdjustTimeStep(adjust);
    }
    if (adjustDelay>-1) {
      auto v = dynamic_cast<VelocityVerletIntegrator*>(integrator);
      if (v) v->setAdjustUpdateDelay(adjustDelay);
    }

    // Check parsing for if any illegal flags were used
    if (check) checkParsing();
  }

  void SimulationBase::standardWriting(bool message) {    
    // Execute writing and saving
#if USE_MPI == 1
    if (rank==0) {
#endif
      if (dataRecord) {
	// Print out time and ratio
	double runTime = dataRecord->getElapsedTime();
	if (!quiet && message) {
	  cout << "Run time: " << runTime << endl;
	  cout << "Ratio: " << dataRecord->getRatio() << endl;
	}

	if (!nowrite) {
	  // Write animation data to files
	  dataRecord->writeData();
	  // Write sectorization data to file
	  dataRecord->writeRunSummary(integrator);
	}

	// Write out stat function data - for now
	if (print) {
	  int numStatFuncs = dataRecord->getNumberOfStatFunctions();
	  for (int i=0; i<numStatFuncs; ++i) {
	    auto data   = dataRecord->getStatFunctionData(i);
	    string name = dataRecord->getStatFunctionName(i);
	    cout << name << "=" << mmPreproc(data) << ";\n";
	    cout << "ListLinePlot[" << name << ",ImageSize->Large,PlotStyle->Black]\n";
	  }
	}
      }
#if USE_MPI == 1
    }
#endif

    // Save arrangement if requested
#if USE_MPI == 1
    if (rank==0) {
#endif
      if (saveFile!="") {
	fileParser.saveToFile(simData, saveFile);
      }
#if USE_MPI == 1
    }
#endif
  }

  void SimulationBase::checkParsing() {
    // Make sure we didn't enter any illegal tokens
    try {
      parser.check();
    }
    catch (ArgParse::UncheckedToken illegal) {
      cout << "Illegal option: [" << illegal.token << "]. Exiting.\n";
      exit(1);
    }
  }

}
