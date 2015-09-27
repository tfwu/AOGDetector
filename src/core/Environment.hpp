/// This file is adapted from UCLA ImageParser by Brandon Rothrock (rothrock@cs.ucla.edu)

#ifndef RGM_ENVIRONMENT_HPP_
#define RGM_ENVIRONMENT_HPP_

#include <vector>

#include "Timer.hpp"
#include "UtilLog.hpp"


namespace RGM
{
/// Predeclaration
class XMLData;

/// Environment for the application
class Environment
{
public:

    Environment();
    ~Environment();

    bool InitFromConfigFile(std::string & configFile);
    bool InitFromConfigString(char * str);

    // Resources
    const std::string&      getConfigFileName() const;
    XMLData*				GetConfig();
    Timers&					GetTimers();

private:
    std::string _configFileName;

    XMLData*    _pConfig;
    Timers*     _pTimers;

    DEFINE_RGM_LOGGER;

}; //class Environment

} // namespace RGM

#endif // RGM_ENVIRONMENT_HPP_


