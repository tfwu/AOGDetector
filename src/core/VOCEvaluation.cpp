#include <fstream>
#include <boost/serialization/vector.hpp>

#include "VOCEvaluation.hpp"
#include "UtilString.hpp"
#include "UtilFile.hpp"
#include "Inference.hpp"


namespace RGM
{

VOCEvaluation::VOCEvaluation(std::string & configFile, int testType)
{
    env_.InitFromConfigFile(configFile);
    XMLData* config = env_.GetConfig();

    switch (testType) {
    case 1:
    {
        detect(config);
        break;
    }
    case 2:
    {
        detectBatch(config);
        break;
    }
    case 0:
    default:
    {
        voc_.read(config);
        data_.setVOCConfig(voc_);
        evalDP();
        break;
    }
    }
}

void VOCEvaluation::evalDP()
{
    const std::vector<std::string>& objs = voc_.ObjectCategories_;

    Inference::Param inferenceParam;

    // original cache dir
    std::string cacheDir = voc_.cacheDir_;

    for ( int i=0; i<objs.size(); ++i ) {
        // Current obj. cat.
        voc_.curObjectIndex_ = i;
        const std::string & objName = voc_.ObjectCategory();

        voc_.cacheDir_ = cacheDir + objName + FileUtil::FILESEP;
        FileUtil::VerifyDirectoryExists(voc_.cacheDir_);

        // check if it is already done
        std::string resultName = voc_.cacheDir_ + "comp3_det_val_" + objName + ".txt";
        if (FileUtil::exists(resultName)) {
            RGM_LOG(normal, boost::format("Testing: %s is already done\n")
                    % objName);
            continue;
        }

        std::string resultNameDet = voc_.cacheDir_ + "comp3_det_val_" + objName + "_dets.bin";
        if ( loadDetections(resultNameDet) ) {
            if ( !data_.prepareTestData(objName) ) {
                continue;
            }
            saveDetBboxes(resultName, data_.testImgFiles());
            continue;
        }

        // check model
        std::string modelName;
        if ( voc_.partConfig_ == VOCConfig::GREEDY_PURSUIT ) {
            //modelName = std::string(objName+"_DPM_final.bin");
            modelName = std::string(objName+"_final.mat.bin");
        } else {
            modelName = std::string(objName+"_RGM_final.bin");
        }

        // Check if it is trained already
        std::string modelFile = voc_.cacheDir_ + modelName;
        if ( !grammar_.read(modelFile) ) {
            RGM_LOG(warning, ("Failed to read model: " + modelFile) );
            continue;
        }
        bool hasPred = grammar_.bboxPred().size();

        // check testing data
        if ( !data_.prepareTestData(objName) ) {
            continue;
        }

        std::string resultNamePt  = voc_.cacheDir_ + "comp3_det_val_" + objName + "_pts.bin";
        if (!loadParseTrees(resultNamePt)) {
            const std::vector<std::string>  & imgFiles = data_.testImgFiles();
            int numImgs = imgFiles.size();
            if (numImgs==0) {
                continue;
            }

            int maxHt = (data_.maxHeight() + grammar_.minCellSize() - 1) / grammar_.minCellSize() + grammar_.pady();
            int maxWd = (data_.maxWidth()  + grammar_.minCellSize() - 1) / grammar_.minCellSize() + grammar_.padx();

            if (!Patchwork::InitFFTW((maxHt + 15) & ~15, (maxWd + 15) & ~15)) {
                RGM_LOG(error, "Error: Could not initialize the Patchwork class\n");
                return;
            }

            // evaluation setting
            grammar_.getInterval() = voc_.interval_;
            Scalar thresh = std::min(voc_.maxThreshold_, grammar_.thresh());
            inferenceParam.useNMS_ = true;
            inferenceParam.nmsOverlap_ = 0.5F;
            inferenceParam.nmsDividedByUnion_ = false;

            // cache FFT
            grammar_.cachingFFTFilters();

            allpts_.clear();
            allpts_.resize(numImgs);

            allDets_.clear();
            allDets_.resize(numImgs);

            Scalar detectLimit = std::numeric_limits<Scalar>::infinity();

#pragma omp parallel for
            for (int j = 0; j < numImgs; ++j ) {
                cv::Mat img;
#pragma omp critical
                {
                    img = cv::imread(imgFiles[j], 1);
                }

                FeaturePyramid pyr(img, grammar_.cellSize(), grammar_.padx(), grammar_.pady(),
                                   0, grammar_.interval(), grammar_.extraOctave());
                if (pyr.empty()) {
                    continue;
                }

                Inference inference(grammar_, inferenceParam);

                if (grammar_.isSingleObjModel()) {
                    inference.runDetection(thresh, pyr, detectLimit, allpts_[j]);
                    for ( int k = 0; k < allpts_[j].size(); ++k ) {
                        if ( hasPred ) {
                            allpts_[j][k].doBboxPred(allDets_[j], k);
                        } else {
                            allpts_[j][k].getSingleObjDet(allDets_[j], k);
                        }
                    }
                } else {
                    inference.runDetectionExt(thresh, pyr, detectLimit, allpts_[j], allDets_[j]);
                    if ( hasPred ) {
                        allDets_[j].clear();
                        for ( int k = 0; k < allpts_[j].size(); ++k ) {
                            allpts_[j][k].doBboxPred(allDets_[j], k);
                        }
                    }
                }

                //#pragma omp critical
                {
                    RGM_LOG(normal, boost::format("  testing %s: %d (%d) ... done.\n")
                            % objName % j % numImgs);
                }
            } // for j

            // save
            saveParseTrees(resultNamePt);
            saveDetections(resultNameDet);

            // write detection resutls for evaluationg AP using PASCAL VOCdevit
            saveDetBboxes(resultName, data_.testImgFiles());

        } else {
            allDets_.clear();
            allDets_.resize(allpts_.size());
            for ( int j = 0; j < allpts_.size(); ++j ) {
                for ( int k = 0; k < allpts_[j].size(); ++k ) {
                    allpts_[j][k].setGrammar(grammar_);
                    if ( hasPred ) {
                        allpts_[j][k].doBboxPred(allDets_[j], k);
                    } else {
                        allpts_[j][k].getSingleObjDet(allDets_[j], k);
                    }
                }
            }
            saveDetections(resultNameDet);           
            // write detection resutls for evaluationg AP using PASCAL VOCdevit
            saveDetBboxes(resultName, data_.testImgFiles());
        }
    } // for i
}

void VOCEvaluation::detect(XMLData* config)
{
    if ( config == NULL ) {
        RGM_LOG(error, "No configuration data specified");
        return;
    }

    std::string modelFile = config->GetString("TestImage/ModelFile");
    float thresh          = config->GetFloat("TestImage/Threshold");
    int interval          = config->GetInt("TestImage/Interval");
    std::string imgFile   = config->GetString("TestImage/ImageFile");

    cv::Mat img = cv::imread(imgFile, cv::IMREAD_COLOR);
    if ( img.empty() ) {
        RGM_LOG(error, boost::format("Can not open image %s") % imgFile);
        return;
    }

    RGM_CHECK(grammar_.read(modelFile), error);

    // Init FFTW
    int maxHt = (img.rows + grammar_.minCellSize() - 1) / grammar_.minCellSize() + grammar_.pady();
    int maxWd = (img.cols + grammar_.minCellSize() - 1) / grammar_.minCellSize() + grammar_.padx();

    if (!Patchwork::InitFFTW((maxHt + 15) & ~15, (maxWd + 15) & ~15)) {
        RGM_LOG(error, "Could not initialize the Patchwork class." );
        return;
    }

    grammar_.getInterval() = interval;
    grammar_.getCachedFFTStatus() = false;

    Inference::Param inferenceParam;
    inferenceParam.useNMS_ = true;
    inferenceParam.nmsOverlap_ = 0.5F;
    inferenceParam.nmsDividedByUnion_ = false;

    //thresh = std::min(grammar_.thresh(), thresh);
    Scalar maxNum = 30000;
    std::vector<ParseTree> pts;

    Inference inference(grammar_, inferenceParam);

    if ( !grammar_.isSingleObjModel() ) {
        std::vector<Detection> dets;
        inference.runDetectionExt(thresh, img, maxNum, pts, dets);

        // bbox prediction
        if ( grammar_.bboxPred().size() > 0 ) {
            std::vector<Detection> preddets;
            for (int i = 0; i < pts.size(); ++i ) {
                pts[i].doBboxPred(preddets, i);
            }
            for ( int i = 0; i < preddets.size(); ++i ) {
                preddets[i].show(img, true);
            }
        } else {
            for ( int i = 0; i < dets.size(); ++i ) {
                dets[i].show(img, true);
            }
        }

    } else {
        inference.runDetection(thresh, img, maxNum, pts);

        // bbox prediction
        if ( grammar_.bboxPred().size() > 0 ) {
            std::vector<Detection> preddets;
            for (int i = 0; i < pts.size(); ++i ) {
                pts[i].doBboxPred(preddets, i);
            }
            for ( int i = 0; i < preddets.size(); ++i ) {
                preddets[i].show(img, true);
            }
        } else {
            for ( int i = 0; i < pts.size(); ++i ) {
                pts[i].showDetection(img, true);
            }
        }
    }
}

void VOCEvaluation::detectBatch(XMLData* config)
{
    if ( config == NULL ) {
        RGM_LOG(error, "No configuration data specified");
        return;
    }

    std::string modelFile = config->GetString("TestBatchImage/ModelFile");
    float thresh          = config->GetFloat("TestBatchImage/Threshold");
    int interval          = config->GetInt("TestBatchImage/Interval");
    std::string imgDir    = config->GetString("TestBatchImage/ImageDir");
    std::string imgExt    = config->GetString("TestBatchImage/ImageExtName");
    int maxImgWd          = config->GetInt("TestBatchImage/maxImageWidth");
    int maxImgHt          = config->GetInt("TestBatchImage/maxImageHeight");
    int numThreads        = config->GetInt("TestBatchImage/numThread");
    numThreads = std::max<int>(numThreads, 1);

    FileUtil::VerifyTheLastFileSep(imgDir);

    // load grammar model
    RGM_CHECK(grammar_.read(modelFile), error);
    std::string objName = grammar_.name();
    if ( objName.empty() ) {
        objName = "NOT_Specified_Obj_Name";
    }

    bool hasPred = grammar_.bboxPred().size();

    // check if it is already done
    std::string resultName = imgDir + "comp3_det_val_" + objName + ".txt";
    if (FileUtil::exists(resultName)) {
        RGM_LOG(normal, boost::format("Testing: %s is already done\n")
                % objName);
        return;
    }

    // get image files
    std::vector<std::string> imgFiles;
    FileUtil::GetFileList(imgFiles, imgDir, imgExt, true);
    if ( imgFiles.size() == 0 ) {
        RGM_LOG(warning, boost::format("Can not find %s images in %s")
                % imgExt % imgDir );
        return;
    }
    int numImg = imgFiles.size();

    std::string resultNameDet = imgDir+ "comp3_det_val_" + objName + "_dets.bin";
    if ( loadDetections(resultNameDet) ) {
        saveDetBboxes(resultName, imgFiles);
        return;
    }

    std::string resultNamePt  = imgDir + "comp3_det_val_" + objName + "_pts.bin";
    if ( loadParseTrees(resultNamePt) ) {
        allDets_.clear();
        allDets_.resize(allpts_.size());
        for ( int j = 0; j < allpts_.size(); ++j ) {
            for ( int k = 0; k < allpts_[j].size(); ++k ) {                
                allpts_[j][k].setGrammar(grammar_);
                if ( hasPred ) {
                    allpts_[j][k].doBboxPred(allDets_[j], k);
                } else {
                    allpts_[j][k].getSingleObjDet(allDets_[j], k);
                }
            }
            std::sort(allDets_[j].begin(), allDets_[j].end());
        }
        saveDetections(resultNameDet);
        // write detection resutls for evaluationg AP using PASCAL VOCdevit
        saveDetBboxes(resultName, imgFiles);
        return;
    }

    RGM_LOG(normal, boost::format("Run batch detection in %s (numImg=%d)") % imgDir % numImg);

    // get the maximum of the sizes of testing images
    if ( maxImgWd <= 0 || maxImgHt <= 0 ) {
       maxImgWd = 0;
       maxImgHt = 0;
       cv::Mat img;
       for ( int i = 0; i < numImg; ++i ) {
           img = cv::imread(imgFiles[i]);
           maxImgWd = std::max<int>(maxImgWd, img.cols);
           maxImgHt = std::max<int>(maxImgHt, img.rows);
       }
    }

    // Init FFTW
    int maxHt = (maxImgHt + grammar_.minCellSize() - 1) / grammar_.minCellSize() + grammar_.pady();
    int maxWd = (maxImgWd + grammar_.minCellSize() - 1) / grammar_.minCellSize() + grammar_.padx();

    if (!Patchwork::InitFFTW((maxHt + 15) & ~15, (maxWd + 15) & ~15)) {
        RGM_LOG(error, "Could not initialize the Patchwork class." );
        return;
    }

    grammar_.getInterval() = interval;
    grammar_.getCachedFFTStatus() = false;

    Inference::Param inferenceParam;
    inferenceParam.useNMS_ = true;
    inferenceParam.nmsOverlap_ = 0.5F;
    inferenceParam.nmsDividedByUnion_ = false;

    thresh = std::min(grammar_.thresh(), thresh);
    Scalar maxNum = 30000; //std::numeric_limits<Scalar>::infinity();;

    allpts_.clear();;
    allpts_.resize(numImg);

    allDets_.clear();
    allDets_.resize(numImg);

    omp_set_dynamic(0);
    // change the number of threads according to the memory and the sizes of images
#pragma omp parallel for num_threads(numThreads)
    for ( int i = 0; i < numImg; ++i ) {
        cv::Mat img;
#pragma omp critical
        {
            img = cv::imread(imgFiles[i], cv::IMREAD_COLOR);
        }

        Inference inference(grammar_, inferenceParam);

        if (grammar_.isSingleObjModel()) {
            inference.runDetection(thresh, img, maxNum, allpts_[i]);
            for ( int k = 0; k < allpts_[i].size(); ++k ) {
                if ( hasPred ) {
                    allpts_[i][k].doBboxPred(allDets_[i], k);
                } else {
                    allpts_[i][k].getSingleObjDet(allDets_[i], k);
                }
            }
        } else {
            inference.runDetectionExt(thresh, img, maxNum, allpts_[i], allDets_[i]);
            if ( hasPred ) {
                allDets_[i].clear();
                for ( int k = 0; k < allpts_[i].size(); ++k ) {
                    allpts_[i][k].doBboxPred(allDets_[i], k);
                }
            }
        }

        RGM_LOG(normal, boost::format("    %d/%d") % i % numImg);
    }

    // save
    saveParseTrees(resultNamePt);
    saveDetections(resultNameDet);

    // write detection resutls for evaluationg AP using PASCAL VOCdevit
    saveDetBboxes(resultName, imgFiles);

    RGM_LOG(normal, "Done.");
}

void VOCEvaluation::saveParseTrees(const std::string & fileName)
{
    std::ofstream out;
    out.open(fileName.c_str(), std::ios::out);
    if ( !out.is_open() ) {
        std::cerr << "Failed to write to file " << fileName << std::endl;
        return;
    }

    boost::archive::binary_oarchive oa(out);

    oa << allpts_;

    out.close();
}

bool VOCEvaluation::loadParseTrees(const std::string & fileName)
{
    std::ifstream in;
    in.open(fileName.c_str(), std::ios::in);
    if ( !in.is_open() ) {
        std::cerr << "Failed to read file " << fileName << std::endl;
        return false;
    }

    allpts_.clear();

    boost::archive::binary_iarchive ia(in);

    ia >> allpts_;

    in.close();

    return true;
}

void VOCEvaluation::saveDetections(const std::string & fileName)
{
    std::ofstream out;
    out.open(fileName.c_str(), std::ios::out);
    if ( !out.is_open() ) {
        std::cerr << "Failed to write to file " << fileName << std::endl;
        return;
    }

    boost::archive::binary_oarchive oa(out);

    oa << allDets_;

    out.close();
}

bool VOCEvaluation::loadDetections(const std::string & fileName)
{
    std::ifstream in;
    in.open(fileName.c_str(), std::ios::in);
    if ( !in.is_open() ) {
        std::cerr << "Failed to read file " << fileName << std::endl;
        return false;
    }

    allDets_.clear();

    boost::archive::binary_iarchive ia(in);

    ia >> allDets_;

    in.close();

    return true;
}

void VOCEvaluation::saveDetBboxes(const std::string &fileName, const std::vector<std::string> &imgFiles)
{
    FILE * f = fopen(fileName.c_str(), "w");
    if (f == NULL) {
        RGM_LOG(error, boost::format("can not write result to %s") % fileName);
        return;
    }

    for ( int j = 0; j < allDets_.size();  ++j ) {
        std::string imageBaseName = FileUtil::GetFileBaseName(imgFiles[j]);
        std::vector<Detection> & det(allDets_[j]);
        //cv::Mat img = cv::imread(imgFiles[j], cv::IMREAD_COLOR);
        for ( int k = 0; k < det.size(); ++k ) {
            const Rectangle_<Scalar>  & r = static_cast<Rectangle_<Scalar> >(det[k]);
            fprintf(f, "%s %f %f %f %f %f\n", // imageName, score, x1, y1, x2, y2
                    imageBaseName.c_str(),
                    det[k].score_, r.left(), r.top(), r.right(), r.bottom());
            //det[k].show(img);
        } // for k
    } // for j

    fclose(f);
}

} // namespace RGM

