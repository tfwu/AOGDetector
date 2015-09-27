#include <sstream>
#include <iomanip> // std::setw

#include "UtilString.hpp"


namespace RGM
{

template<class T>
std::string NumToString_(T t, int zeroPaddingLen)
{
    std::ostringstream  stream;
    if ( zeroPaddingLen <=  0 ) {
        stream << t;
    } else {
        stream << std::setw( zeroPaddingLen ) << std::setfill( '0' ) << t;
    }

    return stream.str();
}

template std::string NumToString_<int>(int t, int zeroPaddingLen);
template std::string NumToString_<float>(float t, int zeroPaddingLen);
template std::string NumToString_<double>(double t, int zeroPaddingLen);


void TokenizeString(std::string str, char delimeter, std::vector<std::string>& tokens)
{
    int offset = 0;
    int pos = str.find(delimeter, offset);
    while (pos != std::string::npos) {
        tokens.push_back(str.substr(offset, pos-offset));
        offset = pos+1;

        pos = str.find(delimeter, offset);
    }
    if (!str.empty() && (offset < str.length()-1)) {
        tokens.push_back(str.substr(offset, str.length()-offset));
    }
}

} // namespace RGM
