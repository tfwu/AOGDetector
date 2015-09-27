/// This file is adapted from UCLA ImageParser by Brandon Rothrock (rothrock@cs.ucla.edu)

#ifndef RGM_COMMANDLINE_HPP_
#define RGM_COMMANDLINE_HPP_

#include <string>
#include <vector>

#include "UtilLog.hpp"

namespace RGM
{
/// Predeclaration
class CommandLineInterpreter;



/// Virtual class for defining commands
class ICommand
{
public:
    virtual std::string GetCommand() = 0;
    virtual int         GetArgumentCount() = 0;
    virtual std::string GetArgumentName(int index) = 0;
    virtual bool        Run(CommandLineInterpreter* pInterpreter) = 0;
protected:
    ICommand() {}

}; //class ICommand


/// Interpret commands
class CommandLineInterpreter
{
public:
    CommandLineInterpreter(int argc, char** argv);

    void AddCommand(ICommand* pCommand);
    int  Process();

    int         GetInt(int index);
    float       GetFloat(int index);
    std::string GetString(int index);

private:
    void PrintHelp();

    std::vector<ICommand*> _commands;
    int                    _argc;
    char**                 _argv;

    DEFINE_RGM_LOGGER;

}; //class CommandLineInterpreter



/// Command for evaluating a model on PASCAL VOC
class VOCEvalCommand : public ICommand
{
public:
    std::string GetCommand();

    int GetArgumentCount();

    std::string GetArgumentName(int index);

    bool Run(CommandLineInterpreter* pInterpreter);

}; //class VOCEvalCommand



/// Command for running detection in a single image
class TestImageCommand : public ICommand
{
public:
    std::string GetCommand();

    int GetArgumentCount();

    std::string GetArgumentName(int index);

    bool Run(CommandLineInterpreter* pInterpreter);

}; //class TestImageCommand


/// Command for running detections in a directory
class TestBatchImageCommand : public ICommand
{
public:
    std::string GetCommand();

    int GetArgumentCount();

    std::string GetArgumentName(int index);

    bool Run(CommandLineInterpreter* pInterpreter);

}; //class TestBatchImageCommand


#ifdef RGM_USE_MATLAB
/// Command for converting DPM VOC-Release5 models
class ConvertDPMVOCRel5ModelCommand : public ICommand
{
public:
    std::string GetCommand();

    int GetArgumentCount();

    std::string GetArgumentName(int index);

    bool Run(CommandLineInterpreter* pInterpreter);
};
#endif


} // namespace RGM

#endif // RGM_COMMANDLINE_HPP_

