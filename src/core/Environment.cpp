#include "Environment.hpp"
#include "UtilFile.hpp"
#include "XMLReader.hpp"

namespace RGM
{

// ---------- Environment ----------

Environment::Environment()
{
    _pConfig = new XMLData();
    _pTimers = new Timers();
}

Environment::~Environment()
{
    if (_pConfig) {
        delete _pConfig;
    }
    if (_pTimers) {
        delete _pTimers;
    }
}

bool Environment::InitFromConfigFile(std::string& configFile)
{
    if (!FileUtil::CheckFileExists(configFile)) {
        RGM_LOG(error, "Config file does not exist:" + configFile);
        return false;
    }

    _configFileName = configFile;

    _pConfig->ReadFromFile(configFile, "Configuration");

    return true;
}

bool Environment::InitFromConfigString(char* str)
{
    _pConfig->ReadFromString(str, "Configuration");

    return true;
}

const std::string& Environment::getConfigFileName() const
{
    return _configFileName;
}

XMLData* Environment::GetConfig()
{
    return _pConfig;
}

Timers& Environment::GetTimers()
{
    return *_pTimers;
}


} // namespace RGM


