#ifndef RGM_UTILSTRING_HPP_
#define RGM_UTILSTRING_HPP_

#include <string>
#include <vector>

namespace RGM
{

/// Convert number to string
template<class T>
std::string NumToString_(T t, int zeroPaddingLen=0);

/// Split a string
void TokenizeString(std::string str, char delimeter, std::vector<std::string>& tokens);

} // namespace RGM

#endif // RGM_UTILSTRING_HPP_
