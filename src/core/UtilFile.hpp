#ifndef RGM_FILEUTILITY_HPP_
#define RGM_FILEUTILITY_HPP_

#include <boost/filesystem.hpp>


namespace RGM
{

namespace fs = boost::filesystem;

/// Utility functions for file management
class FileUtil
{
public:
    /// File separator
#if (defined(WIN32)  || defined(_WIN32) || defined(WIN64) || defined(_WIN64))
    static const char FILESEP = '\\';
#else
    static const char FILESEP = '/';
#endif

    /// check if a file or dir exists
    static bool                     exists(const std::string & input);
    /// Create dir from input ("D:\a\b\c.ext" or "D:\a\b\", i.e. "D:\a\b")
    static void                     CreateDir( const std::string& input );

    /// Get (dir, filebasename, extname) from a full file name
    static std::vector<std::string> FileParts(const std::string& fileName);
    static std::string		        GetParentDir(const std::string& filePath);
    static std::string		        GetFileExtension(const std::string& filePath);
    static std::string		        GetFileBaseName(const std::string& filePath);

    /// Check if a file already existed
    static bool				        CheckFileExists(const std::string &fileName);
    /// Check if a dir existed, if not, create it
    static bool				        VerifyDirectoryExists(const std::string& directoryName, bool create=true);
    /// Delete all files under a dir
    static void				        ClearDirectory(const std::string& directoryName);
    /// Get all files under a dir
    static void				        GetFileList(std::vector<std::string>& files,
                                                const std::string& path,
                                                const std::string& ext, bool useRootDir=false);

    /// Get current work directory
    static std::string              GetCurrentWorkDirectory();

    /// Add the file separator to filePath if necessary
    static void                     VerifyTheLastFileSep(std::string &filePath);

}; // class FileUtil

} // namespace RGM


#endif // RGM_FILEUTILITY_HPP_
