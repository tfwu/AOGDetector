#include <boost/foreach.hpp>

#include "UtilFile.hpp"


namespace RGM
{
bool FileUtil::exists(const std::string & input)
{
    return fs::exists(fs::path(input));
}

void FileUtil::CreateDir( const std::string& input )
{

    fs::path p(input.c_str());

    if (!fs::is_directory(p)) {

        fs::path pp = p.parent_path();
        fs::create_directories(pp);
    }
}

std::vector<std::string> FileUtil::FileParts(const std::string& fileName)
{

    fs::path p( fileName.c_str() );

    std::vector<std::string> vstr(3);

    vstr[0] = p.parent_path().string();
    vstr[1] = p.stem().string();
    vstr[2] = p.extension().string();

    return vstr;
}

std::string FileUtil::GetParentDir(const std::string& filePath)
{

    return fs::path(filePath.c_str()).parent_path().string();
}

std::string FileUtil::GetFileExtension(const std::string& filePath)
{

    return fs::path(filePath.c_str()).extension().string();
}

std::string FileUtil::GetFileBaseName(const std::string& filePath)
{

    return fs::path(filePath.c_str()).stem().string();
}

bool FileUtil::CheckFileExists(const std::string& fileName)
{

    return fs::exists(fs::path(fileName.c_str()));
}

bool FileUtil::VerifyDirectoryExists(const std::string& directoryName, bool create)
{

    fs::path p(directoryName.c_str());

    bool exist = fs::exists(p);

    if (!exist && create) {
        fs::create_directories(p);
    }

    return exist;
}

void FileUtil::ClearDirectory(const std::string& directoryName)
{

    fs::path p(directoryName.c_str());

    fs::remove_all( p );
    fs::create_directory( p );
}

void FileUtil::GetFileList(std::vector<std::string>& files, const std::string& path, const std::string& ext, bool useRootDir)
{

    fs::path targetDir(path.c_str());

    fs::directory_iterator it(targetDir), eod;

    bool getAllFiles = (ext.compare("*")==0);

    BOOST_FOREACH(fs::path const &p, std::make_pair(it, eod)) {

        if (fs::is_regular_file(p)) {
            if ( getAllFiles ) {
                useRootDir ? files.push_back( p.string() ) : files.push_back(p.filename().string());
            } else {
                if (ext.compare(p.extension().string())==0) {
                    useRootDir ? files.push_back( p.string() ) : files.push_back(p.filename().string());
                }
            }
        }
    }
}

std::string FileUtil::GetCurrentWorkDirectory()
{
//#if  (defined(WIN32)  || defined(_WIN32) || defined(WIN64) || defined(_WIN64))
//	char buffer[PATH_MAX];
//	GetModuleFileName( NULL, buffer, PATH_MAX );
//	std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
//	return std::string( buffer ).substr( 0, pos);
//#else
//	char szTmp[32];
//	sprintf(szTmp, "/proc/%d/exe", getpid());
//	int bytes = std::min(readlink(szTmp, pBuf, len), len - 1);
//	if(bytes >= 0)
//		pBuf[bytes] = '\0';
//	return bytes;
//#endif

    fs::path p = fs::initial_path();

    return p.string();
}

void FileUtil::VerifyTheLastFileSep(std::string &filePath)
{
    if ( filePath.find_last_of("/\\") != (filePath.length() - 1) ) {
        filePath = filePath + FILESEP;
    }
}

} // namespace RGM


