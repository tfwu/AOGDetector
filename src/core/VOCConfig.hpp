#ifndef RGM_VOCCONFIG_HPP_
#define RGM_VOCCONFIG_HPP_

#include "XMLReader.hpp"

#include "UtilFile.hpp"
#include "Common.hpp"
#include "UtilLog.hpp"

#include "opencv2/core/types_c.h"

namespace RGM
{
class VOCConfig
{
public:
    /// Methods used to assign latent structures
    enum partConfigType {
        GREEDY_PURSUIT=0, AOG_SEARCH
    };

    /// Default constructor
    VOCConfig();

    /// Reads the configruation settings from a xml file
    void read(XMLData* pConfig);

    /// Finds an object class with given name
    bool findObject(const std::string & objName);

    /// Returns the name of current object category
    const std::string ObjectCategory() const;

    /// Returns the full path of the txt file specifying training POSITIVE dataset
    std::string objectTrainSetFg();

    /// Returns the full path of the txt file specifying training NEGATIVE dataset
    std::string objectTrainSetBg();

    /// Returns the full path of the txt file specifying testing dataset
    std::string objectTestSet();

    /// Data and their directories
    std::vector<std::string> ObjectCategories_;
    int                      curObjectIndex_;

    std::string VOCdevkitDir;
    std::string year_;

    std::string dataDir_; // VOCdevkitDir/VOC+year_
    std::string imgSetDir_;
    std::string annotationDir_;
    std::string imgFileDir_;

    std::string trainSetFg_; // e.g., "trainval"
    std::string trainSetBg_; // e.g., "train"
    std::string testSet_;

    std::string project_; // e.g. "RGM-Release1"

    std::string note_; //decription of the model

    std::string outputDir_;
    std::string modelDir_; // outputDir_/project_/year
    std::string cacheDir_; // outputDir_/project_/year/cache

    // Training
    bool             flipPositives_;
    bool             useDifficultPos_;
    Scalar           C_; // regularization parameter, default 0.001
    Scalar           wlssvmM_;
    Scalar           biasFeature_;

    int              cacheExampleLimit_; // 24000
    int              numNegUsedSmall_;
    int              numNegUsedLarge_;
    float            fgOverlap_;
    int              fgInterval_;
    int              bgInterval_;

    int              dataSplit_;
    int              minClusters_;
    int              maxClusters_;

    int              lrSplitMetric_;

    partConfigType   partConfig_;
    int              partCount_;
    int              partWidth_;
    int              partHeigth_;

    //AOGrid::Param    partConfigParam_;

    int              featSbin_;
    bool             featExtraOctave_;

    // Testing
    int             interval_;
    Scalar          maxThreshold_;

    bool            applyPolicy_;

private:
    DEFINE_RGM_LOGGER;

}; // class VOCConfig

} // namespace RGM

#endif // RGM_VOCCONFIG_HPP_
