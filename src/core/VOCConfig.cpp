#include "VOCConfig.hpp"
#include "UtilString.hpp"

namespace RGM
{

// ---------- VOCConfig ----------

VOCConfig::VOCConfig() :
    curObjectIndex_(-1), flipPositives_(true), useDifficultPos_(false),
    C_(0.001), wlssvmM_(0), biasFeature_(10),
    cacheExampleLimit_(24000), numNegUsedSmall_(200), numNegUsedLarge_(2000),
    fgOverlap_(0.7f), fgInterval_(5), bgInterval_(4),
    dataSplit_(0), lrSplitMetric_(0),
    minClusters_(2), maxClusters_(6), partConfig_(GREEDY_PURSUIT), partCount_(8), partWidth_(6), partHeigth_(6),
    interval_(10), maxThreshold_(-1.1)
{
}

void VOCConfig::read(XMLData *pConfig)
{
    RGM_CHECK_NOTNULL(pConfig);

    // Get directories
    VOCdevkitDir = pConfig->GetString("Paths/VOCdevkitDir");
    FileUtil::VerifyTheLastFileSep(VOCdevkitDir);
    if ( !FileUtil::VerifyDirectoryExists(VOCdevkitDir, false) ) {
        RGM_LOG(error,
                boost::format("not exist VOCdevkitDir %s") % VOCdevkitDir);
        exit(0);
    }

    outputDir_ = pConfig->GetString("Paths/OutputDir");
    FileUtil::VerifyTheLastFileSep(outputDir_);
    FileUtil::VerifyDirectoryExists(outputDir_);

    // Get object categories
    std::string objectCategories = pConfig->GetString("PASCALVOC/ObjectCategories");
    TokenizeString(objectCategories, '+', ObjectCategories_);
    if ( ObjectCategories_.empty() ) {
        RGM_LOG(error, "not specify object categories");
        exit(0);
    }

    year_ = pConfig->GetString("PASCALVOC/Year");
    trainSetFg_ = pConfig->GetString("PASCALVOC/TrainSetFg");
    trainSetBg_ = pConfig->GetString("PASCALVOC/TrainSetBg");
    testSet_ = pConfig->GetString("PASCALVOC/TestSet");

    // Locate image and annotation
    dataDir_ = VOCdevkitDir + "VOC" + year_ + FileUtil::FILESEP;
    if ( !FileUtil::VerifyDirectoryExists(dataDir_) ) {
        RGM_LOG(error,
                boost::format("not exist VOCdataDir %s") % dataDir_);
        exit(0);
    }
    imgSetDir_     = dataDir_ + "ImageSets" + FileUtil::FILESEP + "Main" + FileUtil::FILESEP;
    imgFileDir_    = dataDir_ + "JPEGImages" + FileUtil::FILESEP;
    annotationDir_ = dataDir_ + "Annotations" + FileUtil::FILESEP;

    // Create output directories
    project_ = pConfig->GetString("PASCALVOC/Project");

    note_ = pConfig->GetString("PASCALVOC/Note");

    modelDir_ = outputDir_ + project_ + FileUtil::FILESEP + year_ + FileUtil::FILESEP;
    FileUtil::VerifyDirectoryExists(modelDir_);
    cacheDir_ = modelDir_ + "cache" + FileUtil::FILESEP;
    FileUtil::VerifyDirectoryExists(cacheDir_);

    //////////////////////////////////////////////////////////////
    // Training
    flipPositives_ = pConfig->GetBoolean("TrainParam/FlipPostivies");
    useDifficultPos_ = pConfig->GetBoolean("TrainParam/UseDifficultPostivies");

    C_ = static_cast<Scalar>(pConfig->GetFloat("TrainParam/SVMRegularization"));
    wlssvmM_ = static_cast<Scalar>(pConfig->GetFloat("TrainParam/WLSSVMM"));
    biasFeature_ = static_cast<Scalar>(pConfig->GetFloat("TrainParam/BiasFeature"));

    cacheExampleLimit_ = pConfig->GetInt("TrainParam/CacheExampleLimit");
    numNegUsedSmall_ = pConfig->GetInt("TrainParam/NumUsedNegSmall");
    numNegUsedLarge_ = pConfig->GetInt("TrainParam/NumUsedNegLarge");
    fgOverlap_ = pConfig->GetFloat("TrainParam/FgOverlap");
    fgInterval_ = pConfig->GetInt("TrainParam/FgInterval");
    bgInterval_ = pConfig->GetInt("TrainParam/BgInterval");

    dataSplit_ = pConfig->GetInt("TrainParam/DataSplit");
    lrSplitMetric_ = pConfig->GetInt("TrainParam/DataSplitMetric");
    minClusters_ = pConfig->GetInt("TrainParam/MinClusters");
    maxClusters_ = pConfig->GetInt("TrainParam/MaxClusters");

    partConfig_ = static_cast<partConfigType>(pConfig->GetInt("TrainParam/PartConfiguration"));
    partCount_ = pConfig->GetInt("TrainParam/PartCount");
    partWidth_ = pConfig->GetInt("TrainParam/PartWidth");
    partHeigth_ = pConfig->GetInt("TrainParam/PartHeight");

    featSbin_  = pConfig->GetInt("TrainParam/FeatSbin");
    featExtraOctave_ = pConfig->GetBoolean("TrainParam/FeatExtraOctave");

    //////////////////////////////////////////////////////////////
    // Testing
    interval_ = pConfig->GetInt("TestParam/Interval");
    maxThreshold_ = static_cast<Scalar>(pConfig->GetFloat("TestParam/MaxThreshold"));
    applyPolicy_ = pConfig->GetBoolean("TestParam/ApplyPolicy");    
}

bool VOCConfig::findObject(const std::string & objName)
{
    for ( int i = 0; i < ObjectCategories_.size(); ++i ) {
        if ( objName.compare(ObjectCategories_[i]) == 0 ) {
            curObjectIndex_ = i;
            return true;
        }
    }

    curObjectIndex_ = -1;
    return false;
}

const std::string VOCConfig::ObjectCategory() const
{
    assert(curObjectIndex_ >= 0);
    assert(curObjectIndex_ <  ObjectCategories_.size());
    return ObjectCategories_[curObjectIndex_];
}

std::string VOCConfig::objectTrainSetFg()
{
    return imgSetDir_ + ObjectCategory() + "_" + trainSetFg_ + ".txt";
}

std::string VOCConfig::objectTrainSetBg()
{
    return imgSetDir_ + ObjectCategory() + "_" + trainSetBg_ + ".txt";
}

std::string VOCConfig::objectTestSet()
{
    return imgSetDir_ + ObjectCategory() + "_" + testSet_ + ".txt";
}

} // namespace RGM

