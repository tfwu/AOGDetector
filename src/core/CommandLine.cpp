#include "CommandLine.hpp"

#include <stdio.h>
#include <boost/lexical_cast.hpp>

#include "VOCEvaluation.hpp"

namespace RGM
{

const static char* k_date = __DATE__;
const static char* k_time = __TIME__;
const static char* k_appname = "RGM-Release1.0";


// ---------- CommandLineInterpreter -----------

CommandLineInterpreter::CommandLineInterpreter(int argc, char** argv) :
    _argc(argc), _argv(argv)
{
}

void CommandLineInterpreter::AddCommand(ICommand* pCommand)
{
    _commands.push_back(pCommand);
}

int CommandLineInterpreter::Process()
{
    if (_argc < 2) {
        PrintHelp();
        return 1;
    }

    int argCount = _argc-2;
    std::string command = _argv[1];
    ICommand* pCommand = NULL;
    for (int i=0; i<_commands.size(); i++) {
        if ((_commands[i]->GetCommand().compare(command) == 0) &&
                (_commands[i]->GetArgumentCount() == argCount)) {
            pCommand = _commands[i];
            break;
        }
    }

    if (pCommand == NULL) {
        RGM_LOG(warning, "Unknown command: " + command);
        PrintHelp();
        return 1;
    }

    if (pCommand->GetArgumentCount() != argCount) {
        RGM_LOG(warning, "Wrong number of arguments for command: " + command);
        PrintHelp();
        return 1;
    }

    if (pCommand->Run(this)) {
        return 0;
    }

    return 1;
}

void CommandLineInterpreter::PrintHelp()
{
    printf("%s [built %s %s]\n\nCommands:\n", k_appname, k_date, k_time);
    for (int i=0; i<_commands.size(); i++) {
        printf("%s ", _commands[i]->GetCommand().c_str());
        for (int j=0; j<_commands[i]->GetArgumentCount(); j++) {
            printf("<%s> ", _commands[i]->GetArgumentName(j).c_str());
        }
        printf("\n");
    }
}

int CommandLineInterpreter::GetInt(int index)
{
    return boost::lexical_cast<int>(_argv[index+2]);
}

float CommandLineInterpreter::GetFloat(int index)
{
    return boost::lexical_cast<float>(_argv[index+2]);
}

std::string CommandLineInterpreter::GetString(int index)
{
    return _argv[index+2];
}


// ---------- Commands: VOCEvalCommand ----------

std::string VOCEvalCommand::GetCommand()
{
    return "vocEval";
}

int VOCEvalCommand::GetArgumentCount()
{
    return 1;
}

std::string VOCEvalCommand::GetArgumentName(int index)
{
    switch (index) {
    case 0:
        return "configfile";
    }
    return "";
}

bool VOCEvalCommand::Run(CommandLineInterpreter* pInterpreter)
{
    // Get the config. xml file
    std::string configFileName = pInterpreter->GetString(0);

    VOCEvaluation vocEval(configFileName);

    return true;
}



// ---------- Commands: TestImageCommand ----------

std::string TestImageCommand::GetCommand()
{
    return "detect";
}

int TestImageCommand::GetArgumentCount()
{
    return 1;
}

std::string TestImageCommand::GetArgumentName(int index)
{
    switch (index) {
    case 0:
        return "configfile";
    }
    return "";
}


bool TestImageCommand::Run(CommandLineInterpreter* pInterpreter)
{
    // Get the config. xml file
    std::string configFileName = pInterpreter->GetString(0);

    VOCEvaluation vocEval(configFileName, 1);

    return true;
}

// ---------- Commands: TestBatchImageCommand ----------

std::string TestBatchImageCommand::GetCommand()
{
    return "detectBatch";
}

int TestBatchImageCommand::GetArgumentCount()
{
    return 1;
}

std::string TestBatchImageCommand::GetArgumentName(int index)
{
    switch (index) {
    case 0:
        return "configfile";
    }
    return "";
}


bool TestBatchImageCommand::Run(CommandLineInterpreter* pInterpreter)
{
    // Get the config. xml file
    std::string configFileName = pInterpreter->GetString(0);

    VOCEvaluation vocEval(configFileName, 2);

    return true;
}


#ifdef RGM_USE_MATLAB
// ---------- Commands: ConvertDPMVOCRel5ModelCommand ----------

std::string ConvertDPMVOCRel5ModelCommand::GetCommand()
{
    return "cvtVOCRel5Model";
}

int ConvertDPMVOCRel5ModelCommand::GetArgumentCount()
{
    return 1;
}

std::string ConvertDPMVOCRel5ModelCommand::GetArgumentName(int index)
{
    switch (index) {
    case 0:
        return "vocRel5ModelFile";
    }
    return "";
}

bool ConvertDPMVOCRel5ModelCommand::Run(CommandLineInterpreter* pInterpreter)
{
    // Get the model file
    std::string vocRel5ModelFile = pInterpreter->GetString(0);

    AOGrammar g;

    g.cvtVOCRel5Model(vocRel5ModelFile);

    std::string saveName = vocRel5ModelFile + ".bin";
    g.save(saveName);

    // visualize it
    std::string saveDir = FileUtil::GetParentDir(vocRel5ModelFile);
    g.visualize(saveDir);

    return true;
}

#endif


} // namespace RGM
