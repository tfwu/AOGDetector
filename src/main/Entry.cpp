#include "Common.hpp"
#include "CommandLine.hpp"
#include "UtilLog.hpp"

using namespace RGM;

int main(int argc, char** argv)
{
    log_init();

    DEFINE_RGM_LOGGER;
    RGM_LOG(normal, "Welcome to RGM-release 1.0");

    RGM::CommandLineInterpreter cli(argc, argv);

    cli.AddCommand(new RGM::VOCEvalCommand());    
    cli.AddCommand(new RGM::TestImageCommand());
#ifdef RGM_USE_MATLAB
    cli.AddCommand(new RGM::ConvertDPMVOCRel5ModelCommand());
#endif    
    cli.AddCommand(new TestBatchImageCommand());

    int outcode = cli.Process();

    return outcode;
}
