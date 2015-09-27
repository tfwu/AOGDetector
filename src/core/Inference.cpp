#include <map>
#include <opencv2/core/core_c.h>

#include "Inference.hpp"
#include "UtilGeneric.hpp"

namespace RGM
{

// ------- Inference::Param -------

Inference::Param::Param() :
    thresh_(0.0F), useNMS_(false), nmsOverlap_(0.5F), nmsDividedByUnion_(false),
    createSample_(false), useOverlapLoss_(false), createRootSample2x_(false),
    computeTNodeScores_(false)
{

}

// ------- Inference -------

Inference::Inference(AOGrammar & g, Param & p)
{
    grammar_ = &g;
    param_ = &p;
}

void Inference::runDetection(const Scalar thresh, cv::Mat img, Scalar & maxDetNum,
                             std::vector<ParseTree> & pt, bool usePCA)
{
    if (maxDetNum <= 0) {
        return;
    }

    AOGrammar &g(*grammar_);

    FeaturePyramid pyr(img, g.cellSize(), g.padx(), g.pady(), 0, g.interval(), g.extraOctave());
    if ( pyr.empty() ) {
        return;
    }

    runDetection(thresh, pyr, maxDetNum, pt, usePCA);
}

void Inference::runDetection(const Scalar thresh, const FeaturePyramid & pyramid, Scalar & maxDetNum,
                             std::vector<ParseTree> & pt, bool usePCA)
{
    if (maxDetNum <= 0) {
        return;
    }

    if ( !runDP(pyramid, usePCA) ) {
        return;
    }

    runParsing(thresh, pyramid, maxDetNum, pt, usePCA);
}

void Inference::runDetectionExt(const Scalar thresh, cv::Mat img, Scalar & maxDetNum,
                             std::vector<ParseTree> & pt, std::vector<Detection> & alldets, bool usePCA)
{
    if (maxDetNum <= 0) {
        return;
    }

    AOGrammar &g(*grammar_);

    FeaturePyramid pyr(img, g.cellSize(), g.padx(), g.pady(), 0, g.interval(), g.extraOctave());
    if ( pyr.empty() ) {
        return;
    }

    runDetectionExt(thresh, pyr, maxDetNum, pt, alldets, usePCA);
}

void Inference::runDetectionExt(const Scalar thresh, const FeaturePyramid & pyramid, Scalar & maxDetNum,
                             std::vector<ParseTree> & pt, std::vector<Detection> & alldets, bool usePCA)
{
    if (maxDetNum <= 0) {
        return;
    }

    if ( !runDP(pyramid, usePCA) ) {
        return;
    }

    runParsingExt(thresh, pyramid, maxDetNum, pt, alldets, usePCA);
}

void Inference::runDetection(cv::Mat img, ParseTree &pt, ParseTree &pcaPt)
{
    Scalar thresh = -100;
    Scalar maxDetNum = 1;

    AOGrammar &g(*grammar_);

    RGM_CHECK_NOTNULL(g.pcaCoef());

    FeaturePyramid pyr(img, g.cellSize(), g.padx(), g.pady(), 0, g.interval(), g.extraOctave());
    if ( pyr.empty() ) {
        return;
    }

    // detection for the original model
    std::vector<ParseTree> pts;
    runDetection(thresh, pyr, maxDetNum, pts, false);
    if ( pts.size() == 0 ) {
        return;
    }

    pt = pts[0];
    const ParseInfo * ptInfo = pt.rootNode()->parseInfo(&pt);

//    pt.showDetection(img, true);

    // detection for the pca projected model
    // restricted to the root filter being placed exactly
    // where the original model's root filter was placed
    pyr.project(*g.pcaCoef());

    pyr.getValidLevels().assign(pyr.validLevels().size(), false);
    pyr.getValidLevels()[ptInfo->l_] = true;
    pyr.getValidLevels()[ptInfo->l_ - g.interval()] = true;

    runDP(pyr, true);

    // get the non-loss adjusted detection
    int nbObjComp = grammar_->rootNode()->outEdges().size();
    const Scalar Inf = std::numeric_limits<Scalar>::infinity();
    for ( int c = 0; c < nbObjComp; ++c ) {
        Node * objComp = grammar_->getRootNode()->getOutEdges()[c]->getToNode();
        std::vector<Matrix> & s(getScoreMaps(objComp));
        for ( int l = 0; l < pyr.levels().size(); ++l ) {
            if ( pyr.validLevels()[l] ) {
                if ( c != ptInfo->c_ || l != ptInfo->l_ ) {
                    s[l].setConstant(-Inf);
                } else {
                    Scalar tmp = s[l](ptInfo->y_, ptInfo->x_);
                    s[l].setConstant(-Inf);
                    s[l](ptInfo->y_, ptInfo->x_) = tmp;
                }
            }
        }
    }

    computeORNode(grammar_->getRootNode());

    std::vector<ParseTree> pcaPts;
    runParsing(thresh, pyr, maxDetNum, pcaPts, true);
    RGM_CHECK_EQ(pcaPts.size(), 1);

    pcaPt = pcaPts[0];

//    pcaPt.showDetection(img, true);

    const ParseInfo * pcaPtInfo = pcaPt.rootNode()->parseInfo(&pcaPt);

    RGM_CHECK(pcaPtInfo->c_ == ptInfo->c_ && pcaPtInfo->l_ == ptInfo->l_
              && pcaPtInfo->x_ == ptInfo->x_ && pcaPtInfo->y_ == ptInfo->y_,
              error);
}

bool Inference::runDP(const FeaturePyramid & pyramid, bool usePCA)
{
    // compute score maps for T-nodes
    if (usePCA) {
#if RGM_USE_PCA_DIM
        if (!computePCAAlphaProcesses(pyramid)) {
            RGM_LOG(error, "Failed to computer filter responses." );
            return false;
        }
#else
        RGM_LOG(error, "Not defined RGM_USE_PCA_DIM." );
        return false;
#endif
    } else {
        if ( !computeAlphaProcesses(pyramid) ) {
            RGM_LOG(error, "Failed to computer filter responses." );
            return false;
        }
    }

    computeScalePriorFeature(pyramid.nbLevels());

    // Using DFS order
    std::vector<Node *> & nDFS( grammar_->getNodeDFS() );

    for ( int i = 0; i < nDFS.size(); ++i ) {
        Node * curNode = nDFS[i];
        Node::nodeType t = curNode->type();

        switch ( t ) {
        case Node::AND_NODE: {
            if ( !computeANDNode(curNode, pyramid.padx(), pyramid.pady()) ) {
                return false;
            }

            break;
        }
        case Node::OR_NODE: {
            if ( !computeORNode(curNode) ) {
                return false;
            }

            break;
        }
        } // switch t
    } // for i

    return true;
}

void Inference::runParsing(const Scalar thresh, const FeaturePyramid & pyramid,
                           Scalar & maxDetNum, std::vector<ParseTree> & pt, bool usePCA)
{
    // Find scores above the threshold
    std::vector<Detection> cands;
    for ( int level=0; level<pyramid.nbLevels(); ++level ) {
        if ( !pyramid.validLevels()[level] ) {
            continue;
        }

        const Matrix & score( scoreMaps(grammar_->rootNode())[level] );
        const int rows = score.rows();
        const int cols = score.cols();

        for ( int y=0; y<score.rows(); ++y ) {
            for ( int x=0; x<score.cols(); ++x ) {
                const Scalar s = score(y, x);
                if ( s>thresh ) {
                    // Non-maxima suppresion in a 3x3 neighborhood
                    if (((y == 0) || (x == 0) || (s > score(y - 1, x - 1))) &&
                        ((y == 0) || (s > score(y - 1, x))) &&
                        ((y == 0) || (x == cols - 1) || (s > score(y - 1, x + 1))) &&
                        ((x == 0) || (s > score(y, x - 1))) &&
                        ((x == cols - 1) || (s > score(y, x + 1))) &&
                        ((y == rows - 1) || (x == 0) || (s > score(y + 1, x - 1))) &&
                        ((y == rows - 1) || (s > score(y + 1, x))) &&
                        ((y == rows - 1) || (x == cols - 1) || (s > score(y + 1, x + 1))))
                    {

                        cands.push_back( Detection(level, x, y, s) );	// here, (x, y) is in the coordinate with padding
                    }
                }
            }
        }
    }

    if (cands.empty()) {
        //RGM_LOG(normal, "Not found detections");
        return;
    }

    // Sort scores in descending order
    std::sort(cands.begin(), cands.end());

    if ( cands.size()>maxDetNum ) {
        cands.resize(maxDetNum);
    }

    bool getLoss = lossMaps(grammar_->rootNode()->outEdges()[0]->toNode()).size();

    // Compute detection windows, filter bounding boxes, and derivation trees
    int numDet = cands.size();
    pt.resize(numDet);
    Scalar mInf = -std::numeric_limits<Scalar>::infinity();

#pragma omp parallel for
    for ( int i = 0; i < numDet; ++i ) {
        parse(pyramid, cands[i], pt[i], getLoss, usePCA);
        if ( param_->useNMS_ ) {
            ParseInfo * info = pt[i].getRootNode()->getParseInfo(&pt[i]);
            if ( info->clipBbox(pyramid.imgWd(), pyramid.imgHt()) ) {
                cands[i].c_ = i; // use c_ to record the index which will be used to select the pt after NMS
                cands[i].setX(info->x());
                cands[i].setY(info->y());
                cands[i].setWidth(info->width());
                cands[i].setHeight(info->height());
            } else {
                cands[i].c_ = -1;
                cands[i].score_ = mInf;
            }
        }
    }

    if ( param_->useNMS_ ) {
        std::sort(cands.begin(), cands.end());

        for (int i = 1; i < cands.size(); ++i) {
            cands.resize( std::remove_if(cands.begin() + i, cands.end(),
                                         Intersector_<Scalar>(cands[i - 1], param_->nmsOverlap_, param_->nmsDividedByUnion_)) -
                    cands.begin());
        }

        std::vector<ParseTree> ptNMS;
        ptNMS.reserve(cands.size());

        for ( int i = 0; i < cands.size(); ++i ) {
            if ( cands[i].c_ == -1 ) {
                break;
            }

            int idx = cands[i].c_;
            ptNMS.push_back( pt[idx] );
        }

        pt.swap(ptNMS);
    }
}

void Inference::runParsingExt(const Scalar thresh, const FeaturePyramid & pyramid,
                              Scalar & maxDetNum, std::vector<ParseTree> & pt,
                              std::vector<Detection> & alldets,  bool usePCA)
{
    // Find scores above the threshold
    std::vector<Detection> cands;
    for ( int level=0; level<pyramid.nbLevels(); ++level ) {
        if ( !pyramid.validLevels()[level] ) {
            continue;
        }

        const Matrix & score( scoreMaps(grammar_->rootNode())[level] );
        const int rows = score.rows();
        const int cols = score.cols();

        for ( int y=0; y<score.rows(); ++y ) {
            for ( int x=0; x<score.cols(); ++x ) {
                const Scalar s = score(y, x);
                if ( s>thresh ) {
                    // Non-maxima suppresion in a 3x3 neighborhood
//                    if (((y == 0) || (x == 0) || (s > score(y - 1, x - 1))) &&
//                        ((y == 0) || (s > score(y - 1, x))) &&
//                        ((y == 0) || (x == cols - 1) || (s > score(y - 1, x + 1))) &&
//                        ((x == 0) || (s > score(y, x - 1))) &&
//                        ((x == cols - 1) || (s > score(y, x + 1))) &&
//                        ((y == rows - 1) || (x == 0) || (s > score(y + 1, x - 1))) &&
//                        ((y == rows - 1) || (s > score(y + 1, x))) &&
//                        ((y == rows - 1) || (x == cols - 1) || (s > score(y + 1, x + 1))))
                    {

                        cands.push_back( Detection(level, x, y, s) );	// here, (x, y) is in the coordinate with padding
                    }
                }
            }
        }
    }

    if (cands.empty()) {
        //RGM_LOG(normal, "Not found detections");
        return;
    }

    // Sort scores in descending order
    std::sort(cands.begin(), cands.end());

    if ( cands.size()>maxDetNum ) {
        cands.resize(maxDetNum);
    }

    bool getLoss = lossMaps(grammar_->rootNode()->outEdges()[0]->toNode()).size();

    // Compute detection windows, filter bounding boxes, and derivation trees
    int numDet = cands.size();
    pt.resize(numDet);

#pragma omp parallel for
    for ( int i = 0; i < numDet; ++i ) {
        parse(pyramid, cands[i], pt[i], getLoss, usePCA);
    }

    if ( param_->useNMS_ ) {
        alldets.clear();
        for ( int i = 0; i < pt.size(); ++i ) {
            pt[i].getSingleObjDet(alldets, i);            
        }

        // sort in terms of score
        std::sort(alldets.begin(), alldets.end());

        for ( int i = 1; i < alldets.size(); ++i ) {
            std::vector<Detection> splitdets;

            std::vector<Detection>::iterator iter = alldets.begin() + i;
            for ( ; iter < alldets.end(); ) {
                if ( iter->ptIdx_ == alldets[i-1].ptIdx_ ) {
                    splitdets.push_back(*iter);
                    iter = alldets.erase(iter);
                } else {
                    ++iter;
                }
            }

            alldets.resize( std::remove_if(alldets.begin()+i, alldets.end(),
                                           Intersector_<Scalar>(alldets[i-1],
                                           param_->nmsOverlap_,
                                           param_->nmsDividedByUnion_)) -
                             alldets.begin() );

            alldets.insert(alldets.end(), splitdets.begin(), splitdets.end());

            std::sort(alldets.begin(), alldets.end());
            splitdets.clear();            
        }

        // get all ptIdx not suppressed
        std::vector<int> allptIdx;
        for ( int i = 0; i < alldets.size(); ++i ) {
            allptIdx.push_back(alldets[i].ptIdx_);
        }
        uniqueVector_<int>(allptIdx);

        std::vector<ParseTree> ptNMS;
        ptNMS.reserve(allptIdx.size());

        std::map<int,int> idxmap;
        for ( int i = 0; i < allptIdx.size(); ++i ) {
            ptNMS.push_back( pt[allptIdx[i]] );
            idxmap.insert(std::make_pair(allptIdx[i], i));
        }

        // update the ptIdx of alldets
        for ( int i = 0; i < alldets.size(); ++i ) {
            alldets[i].ptIdx_ = idxmap.find(alldets[i].ptIdx_)->second;
        }


        // update the validality of single Object And nodes        
        for ( int i = 0; i < ptNMS.size(); ++i ) {
            std::vector<PtNode *> sobj = ptNMS[i].getSingleObjAndNodes();
            for ( int j = 0; j < sobj.size(); ++j ) {                
                sobj[j]->getIdx()[PtNode::IDX_VALID] = -1;                
            }
        }      

        for ( int i = 0; i < alldets.size(); ++i ) {            
            std::vector<PtNode *> n = ptNMS[alldets[i].ptIdx_].getNode(alldets[i].c_);
            for ( int j = 0; j < n.size(); ++j ) {
                if ( n[j]->idx()[PtNode::IDX_MYSELF] == alldets[i].ptNodeIdx_ )
                    n[j]->getIdx()[PtNode::IDX_VALID] = 1;
            }
        }

        pt.swap(ptNMS);       
    }
}

void Inference::parse(const FeaturePyramid & pyramid, Detection & cand, ParseTree & pt, bool getLoss, bool usePCA)
{
    pt.clear();
    pt.setGrammar( *grammar_ );

    pt.getImgWd() = pyramid.imgWd();
    pt.getImgHt() = pyramid.imgHt();

    // Backtrack solution in BFS
    std::vector<Node *> gBFS;
    gBFS.push_back( grammar_->getRootNode() );

    // Get the parse info for the root node
    // note that cand.(x_, y_) are in the coordinate with padding
    ParseInfo pinfo(-1, cand.l_, cand.x_, cand.y_, 0, 0, 0, cand.score_, Rectangle_<Scalar>());
    int idxInfo = pt.addParseInfo( pinfo );

    // Add the root node to pt
    int t = static_cast<int>(grammar_->rootNode()->type());
    int gNode = grammar_->idxNode( grammar_->rootNode() );
    pt.getIdxRootNode() = pt.addNode(gNode, t);
    pt.getRootNode()->getIdx()[PtNode::IDX_PARSEINFO] = idxInfo;

    // BFS for pt
    std::vector<int> ptBFS;
    ptBFS.push_back( pt.idxRootNode() );

    int head = 0;
    while ( head < gBFS.size() ) {
        Node * curNode = gBFS[head];
        Node::nodeType t = curNode->type();

        switch ( t ) {
        case Node::T_NODE: {
            if ( !parseTNode(head, gBFS, ptBFS, pyramid, pt, usePCA) ) {
                return;
            }

            break;
        }
        case Node::AND_NODE: {
            if ( !parseANDNode(head, gBFS, ptBFS, pyramid, pt) ) {
                return;
            }

            break;
        }
        case Node::OR_NODE: {
            if ( !parseORNode(head, gBFS, ptBFS, pyramid, pt, getLoss) ) {
                return;
            }

            break;
        }
        default: {
            RGM_LOG(error, "Wrong type of nodes." );
            return;
        }
        }; // switch

        head++;

    } // while
}

bool Inference::computeAlphaProcesses(const FeaturePyramid & pyramid)
{
    // Transform the filters if needed
#pragma omp critical
    if ( !grammar_->cachedFFTStatus() ) {
        grammar_->cachingFFTFilters();
    }

    while (!grammar_->cachedFFTStatus()) {
        RGM_LOG(normal, "Waiting for caching the FFT filters" );
    }

    // Create a patchwork
    const Patchwork patchwork(pyramid);

    // Convolve the patchwork with the filters
    int nbFilters = grammar_->cachedFFTFilters().size();
    std::vector<std::vector<Matrix> > filterResponses(nbFilters); // per Appearance per valid Level

    patchwork.convolve(grammar_->cachedFFTFilters(), filterResponses);

    if ( filterResponses.empty() ) {
        RGM_LOG(error, "filter convolution failed.");
        return false;
    }

    int nbLevel = pyramid.nbLevels();
    int nbValidLevel = filterResponses[0].size();
    assert(nbValidLevel == pyramid.nbValidLevels());

    // score maps of root node
    std::vector<Matrix>& rootScoreMaps( getScoreMaps(grammar_->rootNode()) );
    rootScoreMaps.resize(nbLevel);

    // Normalize the sizes of filter response maps per level
    for ( int l=0, ll=0; l<nbLevel; ++l ) {
        if ( pyramid.validLevels()[l] ) {
            int maxHt = 0;
            int maxWd = 0;
            for ( int i=0; i<nbFilters; ++i ) {
                maxHt = std::max<int>(maxHt, filterResponses[i][ll].rows());
                maxWd = std::max<int>(maxWd, filterResponses[i][ll].cols());
            }

            for ( int i=0; i<nbFilters; ++i ) {
                Matrix tmp = Matrix::Constant(maxHt, maxWd, -std::numeric_limits<Scalar>::infinity());
                tmp.block(0, 0, filterResponses[i][ll].rows(), filterResponses[i][ll].cols()) = filterResponses[i][ll];

                filterResponses[i][ll].swap(tmp);
            }

            rootScoreMaps[l] = Matrix::Zero(maxHt, maxWd);

            ++ll;
        } else {
            rootScoreMaps[l] = Matrix::Zero(1, 1);
        }
    }

    // Assign to T-nodes
    for ( int i=0, t=0; i<grammar_->nodeSet().size(); ++i ) {
        if ( grammar_->nodeSet()[i]->type() == Node::T_NODE ) {
            setScoreMaps( grammar_->nodeSet()[i], nbLevel, filterResponses[t], pyramid.validLevels() );
            ++t;
        }
    }

    return true;
}

#if RGM_USE_PCA_DIM

bool Inference::computePCAAlphaProcesses(const FeaturePyramid & pyramid)
{
    RGM_CHECK_NOTNULL(grammar_);

    if ( pyramid.emptyPCA() ) {
        return false;
    }

    AOGrammar& g(*grammar_);

    // Transform the filters if needed
#pragma omp critical
    if ( !g.cachedFFTStatus() ) {
        g.cachingFFTFilters(true);
    }

    while (!g.cachedFFTStatus()) {
        RGM_LOG(normal, "Waiting for caching the FFT filters" );
    }

    // Create a patchwork
    const Patchwork patchwork(pyramid, true);

    // Convolve the patchwork with the filters
    int nbFilters = g.cachedFFTPCAFilters().size();
    std::vector<std::vector<Matrix> > filterResponses(nbFilters); // per Appearance per valid Level

    patchwork.convolve(g.cachedFFTPCAFilters(), filterResponses);

    if ( filterResponses.empty() ) {
        RGM_LOG(error, "filter convolution failed." );
        return false;
    }

    int nbLevel = pyramid.nbLevels();
    int nbValidLevel = filterResponses[0].size();
    RGM_CHECK_EQ(nbValidLevel, pyramid.nbValidLevels());

    // score maps of root node
    std::vector<Matrix>& rootScoreMaps = getScoreMaps(g.rootNode());
    rootScoreMaps.resize(nbLevel);

    // Normalize the sizes of filter response maps per level
    for ( int l=0, ll=0; l<nbLevel; ++l ) {
        if ( pyramid.validLevels()[l] ) {
            int maxHt = 0;
            int maxWd = 0;
            for ( int i=0; i<nbFilters; ++i ) {
                maxHt = std::max<int>(maxHt, filterResponses[i][ll].rows());
                maxWd = std::max<int>(maxWd, filterResponses[i][ll].cols());
            }

            for ( int i=0; i<nbFilters; ++i ) {
                Matrix tmp = Matrix::Constant(maxHt, maxWd, -std::numeric_limits<Scalar>::infinity());
                tmp.block(0, 0, filterResponses[i][ll].rows(), filterResponses[i][ll].cols()) = filterResponses[i][ll];

                filterResponses[i][ll].swap(tmp);
            }

            rootScoreMaps[l] = Matrix::Zero(maxHt, maxWd);

            ++ll;
        } else {
            rootScoreMaps[l] = Matrix::Zero(1, 1);
        }
    }

    // Assign to T-nodes
    for ( int i=0, t=0; i<g.nodeSet().size(); ++i ) {
        if ( g.nodeSet()[i]->type() == Node::T_NODE ) {
            setScoreMaps( g.nodeSet()[i], nbLevel, filterResponses[t], pyramid.validLevels() );
            ++t;
        }
    }

    return true;
}

#endif // RGM_USE_PCA_DIM

void Inference::computeScalePriorFeature(int nbLevels)
{
    Scaleprior::Param tmp;
    scalepriorFeatures_ = Matrix::Zero(tmp.cols(), nbLevels);

    int s = 0;
    int e = std::min<int>(nbLevels, grammar_->interval());
    scalepriorFeatures_.block(0, s, 1, e).fill(1);

    s = e;
    e = std::min<int>(nbLevels, e*2);
    scalepriorFeatures_.block(1, s, 1, e-s).fill(1);

    s = e;
    scalepriorFeatures_.block(2, s, 1, nbLevels-s).fill(1);
}

bool Inference::computeANDNode(Node * node, int padx, int pady)
{
    if ( node == NULL || node->type() != Node::AND_NODE ) {
        RGM_LOG(error, "Need a valid AND-node as input." );
        return false;
    }

    if ( node->outEdges().size() == 1 && node->outEdges()[0]->type() == Edge::TERMINATION ) {
        return true;
    }

    if ( node->outEdges().size() == 1 && node->outEdges()[0]->type() == Edge::DEFORMATION ) {
        // deformation rule -> apply distance transform
        Deformation::Param w = node->deformationParam();

        // init the score maps using those of the toNode
        std::vector<Matrix>& score( getScoreMaps(node) );
        score = scoreMaps(node->outEdges()[0]->toNode());

        int nbLevel = score.size();

        std::vector<MatrixXi >& x(getDeformationX(node));
        std::vector<MatrixXi >& y(getDeformationY(node));

        x.resize(nbLevel);
        y.resize(nbLevel);

        // Temporary data needed by DT, assume score[0] has the largest size
        int rows = score[0].rows();
        int cols = score[0].cols();

#pragma omp parallel for
        for ( int i=0; i<nbLevel; ++i ) {
            // Bounded distance transform with +/- 4 HOG cells (9x9 window)
            DT2D(score[i], w, Deformation::BoundedShiftInDT, x[i], y[i]);
        }

        return true;
    }

    assert(node->outEdges()[0]->type() == Edge::COMPOSITION);

    // composition rule -> shift and sum scores from toNodes
    std::vector<Matrix>& score( getScoreMaps(node) );
    score = scoreMaps(grammar_->rootNode());

    int nbLevels = score.size();

    // prepare score for this rule
    Scalar bias = 0;
    if ( node->offset() != NULL ) {
        bias = node->offset()->w() * grammar_->featureBias();
    }

    Scaleprior::Vector scalePriorScore = Scaleprior::Vector::Zero(nbLevels);
    if ( node->scaleprior() != NULL ) {
        scalePriorScore = node->scaleprior()->w() * scalepriorFeatures_;
    }

    for ( int i=0; i<nbLevels; ++i ) {
        score[i].fill(bias + scalePriorScore(i));
    }

    Scalar Inf = std::numeric_limits<Scalar>::infinity();

    // sum scores from toNodes (with appropriate shift and down sample)
    std::vector<Edge *> & outEdges = node->getOutEdges();
    for ( int i=0; i<outEdges.size(); ++i ) {
        const Node::Anchor & curAnchor = outEdges[i]->toNode()->anchor();
        int ax = curAnchor(0);
        int ay = curAnchor(1);
        int ds = curAnchor(2);

        // step size for down sampling
        int step = std::pow(2.0f, ds);

        // amount of (virtual) padding to hallucinate
        int virtpady = (step-1)*pady;
        int virtpadx = (step-1)*padx;

        // starting points (simulates additional padding at finer scales)
        // @note (ax, ay) are computed without considering padding.
        // So, given a root location (x, y) in the score map (computed with padded feature map)
        // the location of a part will be: (x-padx) * step + ax without considering padding
        // and (x-padx) * step + ax + padx = x + [ax - (step-1)*padx]
        int starty = ay-virtpady;
        int startx = ax-virtpadx;

        // score table to shift and down sample
        const std::vector<Matrix> & s( scoreMaps(outEdges[i]->toNode()) );

        for (int j = 0; j < s.size(); ++j ) {
            int level = j - grammar_->interval() * ds;
            if (level >= 0 ) {
                // ending points
                int endy = std::min<int>(s[level].rows(), starty + step*(score[j].rows()-1));
                int endx = std::min<int>(s[level].cols(), startx + step*(score[j].cols()-1));

                // y sample points
                std::vector<int> iy;
                int oy = 0;
                for ( int yy=starty; yy<endy; yy+=step ) {
                    if ( yy<0 ) {
                        oy++;
                    } else {
                        iy.push_back(yy);
                    }
                }

                // x sample points
                std::vector<int> ix;
                int ox = 0;
                for ( int xx=startx; xx<endx; xx+=step ) {
                    if ( xx<0 ) {
                        ox++;
                    } else {
                        ix.push_back(xx);
                    }
                }

                // sample scores
                Matrix sp(iy.size(), ix.size());
                for ( int yy=0; yy<iy.size(); ++yy ) {
                    for ( int xx=0; xx<ix.size(); ++xx ) {
                        sp(yy, xx) = s[level](iy[yy], ix[xx]);
                    }
                }

                // sum with correct offset
                Matrix stmp = Matrix::Constant(score[j].rows(), score[j].cols(), -Inf );
                stmp.block(oy, ox, sp.rows(), sp.cols()) = sp;
                score[j] += stmp;

            } else {
                score[j].fill( -Inf );
            }
        }
    }

    return true;
}

void Inference::DT2D(Matrix & scoreMap, Deformation::Param & w, int shift, MatrixXi & Ix, MatrixXi & Iy)
{
    Scalar ax = w(0); // dx^2 dx dy^2 dy
    Scalar bx = w(1);
    Scalar ay = w(2);
    Scalar by = w(3);

    Matrix tmpOut = Matrix::Zero(scoreMap.rows(), scoreMap.cols());
    MatrixXi tmpIy = MatrixXi::Zero(scoreMap.rows(), scoreMap.cols());
    Ix = MatrixXi::Zero(scoreMap.rows(), scoreMap.cols());
    Iy = MatrixXi::Zero(scoreMap.rows(), scoreMap.cols());

    for ( int x=0; x < scoreMap.cols(); ++x ) {
        DT1D(scoreMap.col(x).data(),  tmpOut.col(x).data(), tmpIy.col(x).data(), scoreMap.cols(), shift, scoreMap.rows(), ay, by);
    }

    for ( int y=0; y < scoreMap.rows(); ++y) {
        DT1D(tmpOut.row(y).data(), scoreMap.row(y).data(), Ix.row(y).data(), 1, shift, scoreMap.cols(), ax, bx);
    }

    // get argmax
    for (int x = 0; x < scoreMap.cols(); x++) {
        for (int y = 0; y < scoreMap.rows(); y++) {
            Iy(y, x) = tmpIy(y, Ix(y, x));
        }
    }
}

void Inference::DT1D(const Scalar *vals, Scalar *out_vals, int *I, int step, int shift, int n, Scalar a, Scalar b)
{
    for (int i = 0; i < n; i++) {
        Scalar max_val = -std::numeric_limits<Scalar>::infinity();
        int argmax     = 0;
        int first      = std::max<int>(0,   i-shift);
        int last       = std::min<int>(n-1, i+shift);
        for (int j = first; j <= last; j++) {
            Scalar val = vals[j*step] - a*(i-j)*(i-j) - b*(i-j);
            if (val > max_val) {
                max_val = val;
                argmax  = j;
            }
        }
        out_vals[i*step] = max_val;
        I[i*step] = argmax;
    }
}

bool Inference::computeORNode(Node * node)
{
    if ( node == NULL || node->type() != Node::OR_NODE ) {
        RGM_LOG(error, "Need valid OR-node as input." );
        return false;
    }

    if ( node->outEdges().size()==1 ) {
        return true;
    }

    // take pointwise max over scores of toNodes or outEdges
    std::vector<Matrix>& score( getScoreMaps(node) );
    score = scoreMaps(node->outEdges()[0]->toNode());

    for ( int i = 1; i < node->outEdges().size(); ++i ) {
        for ( int j = 0; j < score.size(); ++j ) {
            score[j] = score[j].cwiseMax(scoreMaps(node->outEdges()[i]->toNode())[j]);
        }
    } // for i

    return true;
}


bool Inference::parseORNode(int idx, std::vector<Node *> & gBFS, std::vector<int> & ptBFS,
                            const FeaturePyramid & pyramid, ParseTree & pt, bool getLoss)
{
    Node * gNode = gBFS[idx];
    if ( gNode->type() != Node::OR_NODE ) {
        RGM_LOG(error, "Not an OR-node." );
        return false;
    }

    int fromIdx = ptBFS[idx];
    PtNode * ptNode = pt.getNodeSet()[fromIdx];

    if (ptNode->getIdx()[PtNode::IDX_PARSEINFO] == -1) {
        RGM_LOG(error, "Need parse info. for the current pt node." );
        return false;
    }

    ParseInfo * info = ptNode->getParseInfo(&pt);

    // Finds the best child of the OR-node by score matching
    int idxArgmax = -1;
    std::vector<Edge *> & outEdges( gNode->getOutEdges() );
    for ( int i = 0; i < outEdges.size(); ++i ) {
        int y = info->y_ - FeaturePyramid::VirtualPadding(pyramid.pady(), info->ds_);
        int x = info->x_ - FeaturePyramid::VirtualPadding(pyramid.padx(), info->ds_);
        Scalar s = scoreMaps(outEdges[i]->toNode())[info->l_](y, x);
        if ( info->score_ == s ) {
            idxArgmax = i;
            break;
        }
    } // for i

    if (idxArgmax == -1) {
        RGM_LOG(error, "Failed to find the best child." );
        return false;
    }

    // Get the best switching
    info->c_ = idxArgmax;
    Edge * bestEdge = outEdges[idxArgmax];
    Node * bestChild = bestEdge->getToNode();

    // Add an edge and a node to pt
    int idxG = grammar_->idxNode(bestChild);
    int t = static_cast<int>(bestChild->type());
    int toIdx = pt.addNode(idxG, t);
    PtNode * toNode = pt.getNodeSet()[toIdx];

    idxG = grammar_->idxEdge(bestEdge);
    int edge = pt.addEdge(fromIdx, toIdx, idxG);

    // Add the node to BFS
    gBFS.push_back( bestEdge->getToNode() );
    ptBFS.push_back( toIdx );

    if ( gNode == grammar_->getRootNode() ) {
        const Rectangle2i & detectWind = bestChild->detectWindow();

        // Detection scale
        Scalar scale = static_cast<Scalar>(grammar_->cellSize()) / pyramid.scales()[info->l_];

        // compute and record image coordinates of the detection window
        Scalar x1 = (info->x_ - pyramid.padx() * std::pow<int>(2, info->ds_)) * scale;
        Scalar y1 = (info->y_ - pyramid.pady() * std::pow<int>(2, info->ds_)) * scale;
        Scalar x2 = x1 + detectWind.width() * scale - 1;
        Scalar y2 = y1 + detectWind.height() * scale - 1;

        // update the parse info.
        info->setX(x1);
        info->setY(y1);
        info->setWidth(x2-x1+1);
        info->setHeight(y2-y1+1);

        if ( getLoss ) {
            const Matrix & lossMap = lossMaps(bestChild)[info->l_];
            info->loss_ = lossMap(info->y_, info->x_);
        }

        // get scale prior and offset feature for toNode
        if ( bestChild->scaleprior() != NULL ) {
            Scaleprior::Param w = scalepriorFeatures_.col(info->l_);
            int idxPrior = pt.addScaleprior(w);
            toNode->getIdx()[PtNode::IDX_SCALEPRIOR] = idxPrior;
        }

        int idxBias = pt.AddBias( grammar_->featureBias() );
        toNode->getIdx()[PtNode::IDX_BIAS] = idxBias;
    }

    // pass the parse info. to the best child
    int idxInfo = pt.addParseInfo(*info); //ptNode->idx()[PtNode::IDX_PARSEINFO]; //
    toNode->getIdx()[PtNode::IDX_PARSEINFO] = idxInfo;

    return true;
}

bool Inference::parseANDNode(int idx, std::vector<Node *> & gBFS, std::vector<int> & ptBFS,
                             const FeaturePyramid & pyramid, ParseTree & pt)
{
    Node * gNode = gBFS[idx];
    if ( gNode->type() != Node::AND_NODE ) {
        RGM_LOG(error, "Not an And-node." );
        return false;
    }

    int fromIdx = ptBFS[idx];
    PtNode * ptNode = pt.getNodeSet()[fromIdx];
    if (ptNode->getIdx()[PtNode::IDX_PARSEINFO] == -1) {
        RGM_LOG(error, "Need parse info. for the current pt node." );
        return false;
    }

    ParseInfo * info = ptNode->getParseInfo(&pt);

    std::vector<Edge *> & outEdges( gNode->getOutEdges() );

    if ( outEdges.size() == 1 && outEdges[0]->type() == Edge::TERMINATION ) {

        // Add an edge and a node to pt
        int idxG = grammar_->idxNode(outEdges[0]->getToNode());
        int t = static_cast<int>(outEdges[0]->getToNode()->type());
        int toIdx = pt.addNode(idxG, t);
        PtNode * toNode = pt.getNodeSet()[toIdx];

        idxG = grammar_->idxEdge(outEdges[0]);
        int edge = pt.addEdge(fromIdx, toIdx, idxG);

        int idxInfo = pt.addParseInfo(*info); //ptNode->idx()[PtNode::IDX_PARSEINFO]; //
        toNode->getIdx()[PtNode::IDX_PARSEINFO] = idxInfo;

        // Add the node to BFS
        gBFS.push_back(outEdges[0]->getToNode());
        ptBFS.push_back(toIdx);

        return true;
    }

    if ( outEdges.size() == 1 && outEdges[0]->type() == Edge::DEFORMATION ) {

        const MatrixXi & Ix =  getDeformationX(gNode)[info->l_];
        const MatrixXi & Iy =  getDeformationY(gNode)[info->l_];

        const int vpadx = FeaturePyramid::VirtualPadding(pyramid.padx(), info->ds_);
        const int vpady = FeaturePyramid::VirtualPadding(pyramid.pady(), info->ds_);

        // Location of ptNode without virtual padding
        int nvpX = info->x_ - vpadx;
        int nvpY = info->y_ - vpady;

        // Computing the toNode's location:
        //  - the toNode is (possibly) deformed to some other location
        //  - lookup its displaced location using the distance transform's argmax tables Ix and Iy
        int defX = Ix(nvpY, nvpX);
        int defY = Iy(nvpY, nvpX);

        // with virtual padding
        int toX = defX + vpadx;
        int toY = defY + vpady;

        // get deformation vectors
        int dx = info->x_ - toX;
        int dy = info->y_ - toY;

        if (ptNode->idx()[PtNode::IDX_DEF] != -1) {
            RGM_LOG(error, "Parsing wrong deformation AND-node" );
            return false;
        }

        if (ptNode->idx()[PtNode::IDX_DEF] == -1 ) {
            ptNode->getIdx()[PtNode::IDX_DEF] = pt.addDeformation(dx, dy, gNode->isLRFlip());
        }

        // Look up the score of toNode
        const Matrix & score = scoreMaps(outEdges[0]->toNode())[info->l_];
        Scalar s = score(defY, defX);

        // Add an edge and a node to pt
        int idxG = grammar_->idxNode(outEdges[0]->getToNode());
        int t = static_cast<int>(outEdges[0]->getToNode()->type());
        int toIdx = pt.addNode(idxG, t);
        PtNode * toNode = pt.getNodeSet()[toIdx];

        idxG = grammar_->idxEdge(outEdges[0]);
        int edge = pt.addEdge(fromIdx, toIdx, idxG);

        // Detection scale and window
        Scalar scale = static_cast<Scalar>(grammar_->cellSize()) / pyramid.scales()[info->l_];
        const Rectangle2i & detectWind = outEdges[0]->toNode()->detectWindow();

        // compute and record image coordinates of the detection window
        Scalar x1 = (toX - pyramid.padx() * std::pow<int>(2, info->ds_)) * scale;
        Scalar y1 = (toY - pyramid.pady() * std::pow<int>(2, info->ds_)) * scale;
        Scalar x2 = x1 + detectWind.width()*scale - 1;
        Scalar y2 = y1 + detectWind.height()*scale - 1;

        ParseInfo pinfo(0, info->l_, toX, toY, info->ds_, dx, dy, s,
                        Rectangle_<Scalar>(x1, y1, x2-x1+1, y2-y1+1));

        toNode->getIdx()[PtNode::IDX_PARSEINFO] = pt.addParseInfo( pinfo );

        // Add the node to BFS
        gBFS.push_back(outEdges[0]->getToNode());
        ptBFS.push_back(toIdx);

        return true;
    }

    RGM_CHECK(outEdges.size() >= 1 && outEdges[0]->type() == Edge::COMPOSITION, error);
    for ( int i = 0; i < outEdges.size(); ++i ) {
        // get anchor
        const Node::Anchor & anchor = outEdges[i]->toNode()->anchor();
        int ax = anchor(0);
        int ay = anchor(1);
        int ads = anchor(2);

        // compute the location of toNode
        int toX = info->x_ * std::pow<int>(2, ads) + ax;
        int toY = info->y_ * std::pow<int>(2, ads) + ay;
        int toL = info->l_ - grammar_->interval() * ads;

        // Accumulate rescalings relative to ptNode
        int tods = info->ds_ + ads;

        // get the score of toNode
        const Matrix & score = scoreMaps(outEdges[i]->toNode())[toL];
        int nvpX = toX - FeaturePyramid::VirtualPadding(pyramid.padx(), tods);
        int nvpY = toY - FeaturePyramid::VirtualPadding(pyramid.pady(), tods);
        Scalar s = score(nvpY, nvpX);

        // Detection scale and window
        Scalar scale = static_cast<Scalar>(grammar_->cellSize()) / pyramid.scales()[toL];
        const Rectangle2i & detectWind = outEdges[i]->toNode()->detectWindow();

        // compute and record image coordinates of the detection window
        Scalar x1 = (toX - pyramid.padx()* std::pow<int>(2, tods)) * scale;
        Scalar y1 = (toY - pyramid.pady()* std::pow<int>(2, tods)) * scale;
        Scalar x2 = x1 + detectWind.width()*scale - 1;
        Scalar y2 = y1 + detectWind.height()*scale - 1;

        // Add an edge and a node to pt
        int idxG = grammar_->idxNode(outEdges[i]->getToNode());
        int t = static_cast<int>(outEdges[i]->getToNode()->type());
        int toIdx = pt.addNode(idxG, t);
        PtNode * toNode = pt.getNodeSet()[toIdx];

        idxG = grammar_->idxEdge(outEdges[i]);
        int edge = pt.addEdge(fromIdx, toIdx, idxG);

        ParseInfo pinfo(0, toL, toX, toY, tods, 0, 0, s,
                        Rectangle_<Scalar>(x1, y1, x2-x1+1, y2-y1+1));

        toNode->getIdx()[PtNode::IDX_PARSEINFO] = pt.addParseInfo(pinfo);

        // Add the node to BFS
        gBFS.push_back(outEdges[i]->getToNode());
        ptBFS.push_back(toIdx);

    } // for i

    return true;
}

bool Inference::parseTNode(int idx, std::vector<Node *> & gBFS, std::vector<int> & ptBFS,
                           const FeaturePyramid & pyramid, ParseTree & pt, bool usePCA)
{
    Node * gNode = gBFS[idx];
    if ( gNode->type() != Node::T_NODE ) {
        RGM_LOG(error, "Not an T-node." );
        return false;
    }

    int fromIdx = ptBFS[idx];
    PtNode * ptNode = pt.getNodeSet()[fromIdx];
    if (ptNode->getIdx()[PtNode::IDX_PARSEINFO] == -1) {
        RGM_LOG(error, "Need parse info. for the current pt node." );
        return false;
    }

    if ( param_->createSample_ ) {
        ParseInfo * info = ptNode->getParseInfo(&pt);

        int ds2 = std::pow<int>(2, info->ds_) - 1;
        int fy = info->y_ - pyramid.pady() * ds2;
        int fx = info->x_ - pyramid.padx() * ds2;
        Appearance::Param w = pyramid.levels()[info->l_].block(fy, fx, gNode->appearance()->w().rows(), gNode->appearance()->w().cols());

        //int bs = 20;
        /*cv::Mat img = FeaturePyramid::visualize(pyramid.levels()[info->l_], bs);
        cv::rectangle(img, cv::Rect(fx*bs, fy*bs, gNode->appearance()->w().cols()*bs, gNode->appearance()->w().rows()*bs), cv::Scalar::all(255), 3);
        cv::imshow("HOG", img);
        cv::waitKey(0);

        if ( !getTestImg().empty() ) {
            cv::rectangle(getTestImg(), info->cvRect(), cv::Scalar::all(0), 2);
            cv::imshow("create sample", getTestImg());
            cv::waitKey(0);
        }*/

        //FeaturePyramid::visualize(w, bs);

        ptNode->getIdx()[PtNode::IDX_APP] = pt.addAppearance(w, gNode->isLRFlip());

        if ( pt.appUsage().size() == 0 ) {
            pt.getAppUsage().assign(grammar_->appearanceSet().size(), 0);
        }

        pt.getAppUsage()[grammar_->idxAppearance(gNode->appearance())] += 1;
    }

    if ( param_->createRootSample2x_ ) {
        ParseInfo * info = ptNode->getParseInfo(&pt);

        int ads = 1;
        int ax = 0;
        int ay = 0;

        int ads2 = std::pow<int>(2, ads);

        int toX = info->x_ * ads2 + ax;
        int toY = info->y_ * ads2 + ay;
        int toL = info->l_ - grammar_->interval() * ads;

        int ds2 = ads2 - 1;
        int fy = toY - pyramid.pady() * ds2;
        int fx = toX - pyramid.padx() * ds2;

        int wd = gNode->appearance()->w().rows() * ads2;
        int ht = gNode->appearance()->w().cols() * ads2;

        assert(pt.appearaceX() == NULL);
        pt.getAppearaceX() = new Appearance::Param();

        if ( (fy >= 0) && (fy + ht <= pyramid.levels()[toL].rows()) &&
             (fx >= 0) && (fx + wd <= pyramid.levels()[toL].cols()) ) {
            *pt.getAppearaceX() = pyramid.levels()[toL].block(fy, fx, ht, wd);
        } else {
            *pt.getAppearaceX() = Appearance::Param::Constant( ht, wd, pyramid.levels()[toL](0, 0) );

                    int x1 = std::max<int>(fx, 0);
            int x2 = std::min<int>(fx + wd, pyramid.levels()[toL].cols());
            int y1 = std::max<int>(fy, 0);
            int y2 = std::min<int>(fy + ht, pyramid.levels()[toL].rows());
            int wd1 = x2 - x1;
            int ht1 = y2 - y1;

            int fx2 = (fx >= 0) ? 0 : -fx;
            int fy2 = (fy >= 0) ? 0 : -fy;

            pt.getAppearaceX()->block(fy2, fx2, ht1, wd1) = pyramid.levels()[toL].block(y1, x1, ht1, wd1);
        }
        if (gNode->isLRFlip()) {
            *pt.getAppearaceX() = FeaturePyramid::Flip(*pt.getAppearaceX());
        }

        if ( pt.appUsage().size() == 0 ) {
            pt.getAppUsage().assign(grammar_->appearanceSet().size(), 0);
        }

        pt.getAppUsage()[grammar_->idxAppearance(gNode->appearance())] += 1;

        //int bs = 20;
        //FeaturePyramid::visualize(*pt.appearaceX(), bs);

        //Appearance::Param w = pyramid.levels()[info->l_].block(info->y_, info->x_, gNode->appearance()->w().rows(), gNode->appearance()->w().cols());
        //FeaturePyramid::visualize(w, bs);
    }

    if ( param_->computeTNodeScores_ ) {
        ParseInfo * info = ptNode->getParseInfo(&pt);

        int ds2 = std::pow<int>(2, info->ds_) - 1;
        int fy = info->y_ - pyramid.pady() * ds2;
        int fx = info->x_ - pyramid.padx() * ds2;
        if ( usePCA ) {
#if RGM_USE_PCA_DIM
            const FeaturePyramid::PCALevel feat = pyramid.PCAlevels()[info->l_].block(fy, fx, gNode->wpca()->rows(), gNode->wpca()->cols());

            Eigen::Map<const Matrix, Eigen::Aligned> mapF(feat.data()->data(),
                gNode->wpca()->rows(), gNode->wpca()->cols() * FeaturePyramid::NbPCAFeatures);

            Eigen::Map<const Matrix, Eigen::Aligned> mapW(gNode->getWpca()->data()->data(),
                gNode->wpca()->rows(), gNode->wpca()->cols() * FeaturePyramid::NbPCAFeatures);

            info->score_ = (mapF.cwiseProduct(mapW)).sum();
#endif
        } else {
            const FeaturePyramid::Level feat = pyramid.levels()[info->l_].block(fy, fx, gNode->appearance()->w().rows(), gNode->appearance()->w().cols());

            Eigen::Map<const Matrix, Eigen::Aligned> mapF(feat.data()->data(),
                gNode->appearance()->w().rows(), gNode->appearance()->w().cols() * FeaturePyramid::NbFeatures);

            Appearance::Param w = gNode->appearanceParam();
            Eigen::Map<const Matrix, Eigen::Aligned> mapW(w.data()->data(),
                w.rows(), w.cols() * FeaturePyramid::NbFeatures);

            info->score_ = (mapF.cwiseProduct(mapW)).sum();
        }
    }

    return true;
}

const std::vector<Matrix >& Inference::scoreMaps(const Node * n)
{
    RGM_CHECK_NOTNULL(n);

    if ( n->type() == Node::OR_NODE &&
         n->outEdges().size() == 1 &&
         n->outEdges()[0]->type() == Edge::SWITCHING )  {

        return scoreMaps(n->outEdges()[0]->toNode());
    }

    if (n->type() == Node::AND_NODE &&
            n->outEdges().size() == 1 &&
            n->outEdges()[0]->type() == Edge::TERMINATION ) {

        return scoreMaps(n->outEdges()[0]->toNode());
    }

    Maps::const_iterator iter = scoreMaps_.find(n->tag());
    RGM_CHECK_NOTEQ( iter,  scoreMaps_.end() );

    return iter->second;
}

std::vector<Matrix >& Inference::getScoreMaps(const Node * n)
{
    RGM_CHECK_NOTNULL(n);

    if ( n->type() == Node::OR_NODE &&
         n->outEdges().size() == 1 &&
         n->outEdges()[0]->type() == Edge::SWITCHING )  {

        return getScoreMaps(n->outEdges()[0]->toNode());
    }

    if (n->type() == Node::AND_NODE &&
            n->outEdges().size() == 1 &&
            n->outEdges()[0]->type() == Edge::TERMINATION ) {

        return getScoreMaps(n->outEdges()[0]->toNode());
    }

    Maps::iterator iter = scoreMaps_.find(n->tag());
    if ( iter == scoreMaps_.end() ) {
        scoreMaps_.insert(std::make_pair(n->tag(), std::vector<Matrix >()));
    }

    return scoreMaps_[n->tag()];
}

void Inference::setScoreMaps(const Node * n, int nbLevels, std::vector<Matrix> & s, const std::vector<bool> & validLevels)
{
    RGM_CHECK_NOTNULL(n);

    if (n->type() != Node::T_NODE ) {
        return;
    }

    std::vector<Matrix >& m( getScoreMaps(n) );

    m.resize(nbLevels);

    for ( int i=0, j=0; i<nbLevels; ++i) {
        if ( validLevels[i] ) {
            m[i].swap(s[j]);
            ++j;
        } else {
            m[i] = Matrix::Constant(1, 1, -std::numeric_limits<Scalar>::infinity());
        }
    }
}

const std::vector<Matrix >& Inference::scoreMapCopies(const Node * n)
{
    RGM_CHECK_NOTNULL(n);

    if ( n->type() == Node::OR_NODE &&
         n->outEdges().size() == 1 &&
         n->outEdges()[0]->type() == Edge::SWITCHING )  {

        return scoreMapCopies(n->outEdges()[0]->toNode());
    }

    if (n->type() == Node::AND_NODE &&
            n->outEdges().size() == 1 &&
            n->outEdges()[0]->type() == Edge::TERMINATION ) {

        return scoreMapCopies(n->outEdges()[0]->toNode());
    }

    Maps::const_iterator iter = scoreMapCopies_.find(n->tag());
    RGM_CHECK_NOTEQ( iter, scoreMapCopies_.end() );

    return iter->second;
}

std::vector<Matrix >& Inference::getScoreMapCopies(const Node * n)
{
    RGM_CHECK_NOTNULL(n);

    if ( n->type() == Node::OR_NODE &&
         n->outEdges().size() == 1 &&
         n->outEdges()[0]->type() == Edge::SWITCHING )  {

        return getScoreMapCopies(n->outEdges()[0]->toNode());
    }

    if (n->type() == Node::AND_NODE &&
            n->outEdges().size() == 1 &&
            n->outEdges()[0]->type() == Edge::TERMINATION ) {

        return getScoreMapCopies(n->outEdges()[0]->toNode());
    }

    Maps::iterator iter = scoreMapCopies_.find(n->tag());
    if ( iter == scoreMapCopies_.end() ) {
        scoreMapCopies_.insert(std::make_pair(n->tag(), std::vector<Matrix >()));
    }

    return scoreMapCopies_[n->tag()];
}

const std::vector<bool>& Inference::scoreMapStatus(const Node * n)
{

    //RGM_CHECK_NOTNULL(n);

    if ( n->type() == Node::OR_NODE &&
         n->outEdges().size() == 1 &&
         n->outEdges()[0]->type() == Edge::SWITCHING )  {

        return scoreMapStatus(n->outEdges()[0]->toNode());
    }

    if (n->type() == Node::AND_NODE &&
            n->outEdges().size() == 1 &&
            n->outEdges()[0]->type() == Edge::TERMINATION ) {

        return scoreMapStatus(n->outEdges()[0]->toNode());
    }

    Status::const_iterator iter = scoreMapStatus_.find(n->tag());
    RGM_CHECK_NOTEQ( iter, scoreMapStatus_.end() );

    return iter->second;
}

std::vector<bool>& Inference::getScoreMapStatus(const Node * n)
{
    //RGM_CHECK_NOTNULL(n);

    if ( n->type() == Node::OR_NODE &&
         n->outEdges().size() == 1 &&
         n->outEdges()[0]->type() == Edge::SWITCHING )  {

        return getScoreMapStatus(n->outEdges()[0]->toNode());
    }

    if (n->type() == Node::AND_NODE &&
            n->outEdges().size() == 1 &&
            n->outEdges()[0]->type() == Edge::TERMINATION ) {

        return getScoreMapStatus(n->outEdges()[0]->toNode());
    }

    Status::iterator iter = scoreMapStatus_.find(n->tag());
    if ( iter == scoreMapStatus_.end() ) {
        scoreMapStatus_.insert(std::make_pair(n->tag(), std::vector<bool >()));
    }

    return scoreMapStatus_[n->tag()];
}

void Inference::setScoreMapStatus(const Node* n, int l)
{
    int num = scoreMaps(n).size();

    RGM_CHECK_GE(l, 0);
    RGM_CHECK_LT(l, num);

    std::vector<bool>& status(getScoreMapStatus(n));
    if ( num != status.size() ) {
        status.assign(num, false);
    }

    status[l] = true;

}

void Inference::copyScoreMaps(const Node * n)
{
    const std::vector<Matrix >& s(scoreMaps(n));

    std::vector<Matrix >& sCopies(getScoreMapCopies(n));
    sCopies.resize(s.size());

    std::copy(s.begin(), s.end(), sCopies.begin());

    // set score map status
    std::vector<bool>& status(getScoreMapStatus(n));
    status.assign(s.size(), false);
}

void Inference::recoverScoreMaps(const Node * n)
{
    const std::vector<bool>& status(scoreMapStatus(n));

    for ( int i = 0; i < status.size(); ++i ) {
        if ( status[i] ) {
            getScoreMaps(n)[i] = scoreMapCopies(n)[i];
            getScoreMapStatus(n)[i] = false;
        }
    }
}

std::vector<MatrixXi>& Inference::getDeformationX(const Node * n)
{
    RGM_CHECK_NOTNULL(n);

    argMaps::iterator iter = deformationX_.find(n->tag());
    if ( iter == deformationX_.end() ) {
        deformationX_.insert(std::make_pair(n->tag(), std::vector<MatrixXi >()));
    }

    return deformationX_[n->tag()];
}

std::vector<MatrixXi>& Inference::getDeformationY(const Node * n)
{
    RGM_CHECK_NOTNULL(n);

    argMaps::iterator iter = deformationY_.find(n->tag());
    if ( iter == deformationY_.end() ) {
        deformationY_.insert(std::make_pair(n->tag(), std::vector<MatrixXi >()));
    }

    return deformationY_[n->tag()];
}


const std::vector<Matrix >& Inference::lossMaps(const Node * n)
{
    RGM_CHECK_NOTNULL(n);

    if ( n->type() == Node::OR_NODE &&
         n->outEdges().size() == 1 &&
         n->outEdges()[0]->type() == Edge::SWITCHING )  {

        return lossMaps(n->outEdges()[0]->toNode());
    }

    if (n->type() == Node::AND_NODE &&
            n->outEdges().size() == 1 &&
            n->outEdges()[0]->type() == Edge::TERMINATION ) {

        return lossMaps(n->outEdges()[0]->toNode());
    }

    Maps::const_iterator iter = lossMaps_.find(n->tag());

    if ( iter != lossMaps_.end() ) {
        return iter->second;
    } else {
        return std::vector<Matrix>();
    }

}

std::vector<Matrix >& Inference::getLossMaps(const Node * n)
{
    RGM_CHECK_NOTNULL(n);

    if ( n->type() == Node::OR_NODE &&
         n->outEdges().size() == 1 &&
         n->outEdges()[0]->type() == Edge::SWITCHING )  {

        return getLossMaps(n->outEdges()[0]->toNode());
    }

    if (n->type() == Node::AND_NODE &&
            n->outEdges().size() == 1 &&
            n->outEdges()[0]->type() == Edge::TERMINATION ) {

        return getLossMaps(n->outEdges()[0]->toNode());
    }

    Maps::iterator iter = lossMaps_.find(n->tag());
    if ( iter == lossMaps_.end() ) {
        lossMaps_.insert(std::make_pair(n->tag(), std::vector<Matrix >()));
    }

    return lossMaps_[n->tag()];
}


std::vector<bool> Inference::computeOverlapMaps(const std::vector<Rectangle2i> & bboxes, const FeaturePyramid & pyr,
                                                std::vector<std::vector<std::vector<Matrix> > > & overlapMaps, Scalar overlapThr)
{
    int nbObjComp = grammar_->rootNode()->outEdges().size();

    const std::vector<Matrix> & s(scoreMaps(grammar_->rootNode()));

    overlapMaps.resize(s.size());

    std::vector<bool> valid(s.size(), false);

    for ( int l = 0; l < s.size(); ++l ) {
        overlapMaps[l].resize(bboxes.size());

        const int rows = s[l].rows();
        const int cols = s[l].cols();
        const Scalar scale = static_cast<Scalar>(grammar_->cellSize()) / pyr.scales()[l];

        for ( int b = 0; b < bboxes.size(); ++b ) {
            overlapMaps[l][b].resize(nbObjComp);

            Rectangle_<float> refBox(bboxes[b].x(), bboxes[b].y(), bboxes[b].width(), bboxes[b].height()); // at original image resolution
            Intersector_<float> inter(refBox, overlapThr, true);

            const bool imgClip = (bboxes[b].area() / float(pyr.imgWd()*pyr.imgHt())) < 0.7F;

            for ( int c = 0; c < nbObjComp; ++c ) {
                const Rectangle2i & detWind = grammar_->rootNode()->outEdges()[c]->toNode()->detectWindow();
                const float wd = detWind.width() * scale;
                const float ht = detWind.height() * scale;
                assert(wd > 0 && ht >0);

                Matrix & o = overlapMaps[l][b][c];
                o = Matrix::Zero(rows, cols);

                for ( int y = 0; y < std::max<int>(o.rows(), pyr.levels()[l].rows()); ++y ) {
                    float y1 = (y - grammar_->pady()) * scale;
                    float y2 = y1 + ht - 1;
                    if ( imgClip ) {
                        y1 = std::min<float>(pyr.imgHt()-1, std::max<float>(y1, 0));
                        y2 = std::min<float>(pyr.imgHt()-1, std::max<float>(y2, 0));
                    }

                    for ( int x = 0; x < std::max<int>(o.cols(), pyr.levels()[l].cols()); ++x ) {
                        float x1 = (x - grammar_->padx()) * scale;
                        float x2 = x1 + wd - 1;
                        if ( imgClip ) {
                            x1 = std::min<float>(pyr.imgWd()-1, std::max<float>(x1, 0));
                            x2 = std::min<float>(pyr.imgWd()-1, std::max<float>(x2, 0));
                        }

                        Rectangle_<float> box(x1, y1, x2-x1+1, y2-y1+1);
                        float ov = 0;

                        if ( inter(box, &ov) && (y<pyr.levels()[l].rows()) && (x<pyr.levels()[l].cols())) {
                            valid[l] = true;
                            // assumes that models only have one level of parts
                            if (l-pyr.interval()>=0) {
                                valid[l-pyr.interval()] = true;
                            }
                        }

                        if ( y < o.rows() && x < o.cols() ) {
                            o(y, x) = ov;
                        }

                    } // for x
                } // for y
            } // for c
        } // for b
    } // for l

    for ( int l = 0; l < s.size(); ++l ) {
        if ( !valid[l] ) {
            overlapMaps[l].clear();
        }
    }

    return valid;
}


void Inference::inhibitOutput(int idxBox, std::vector<std::vector<std::vector<Matrix> > > & overlapMap, Scalar overlapThr, bool needCopy)
{
    int nbObjComp = grammar_->rootNode()->outEdges().size();

    for ( int c = 0; c < nbObjComp; ++c ) {
        Node * objComp = grammar_->getRootNode()->getOutEdges()[c]->getToNode();
        std::vector<Matrix> & s(getScoreMaps(objComp));

        for ( int i = 0; i < s.size(); ++i ) {
            if ( overlapMap[i].size() == 0 ) {
                continue;
            }

            if ( needCopy ) {
                setScoreMapStatus(objComp, i);
            }

            Matrix & o = overlapMap[i][idxBox][c];
            s[i] = (o.array() < overlapThr).select(-std::numeric_limits<Scalar>::infinity(), s[i]);
        }
    } // for c

    computeORNode(grammar_->getRootNode());
}


void Inference::applyLossAdjustment(int idxBox, int nbBoxes, std::vector<std::vector<std::vector<Matrix> > > & overlapMap,
                                    Scalar fgOverlap, Scalar bgOverlap, bool needCopy)
{

    int nbObjComp = grammar_->rootNode()->outEdges().size();

    const Scalar Inf = std::numeric_limits<Scalar>::infinity();

    for ( int c = 0; c < nbObjComp; ++c ) {
        Node * objComp = grammar_->getRootNode()->getOutEdges()[c]->getToNode();

        std::vector<Matrix> & s(getScoreMaps(objComp));

        std::vector<Matrix> & lmap(getLossMaps(objComp));
        lmap.resize(s.size(), Matrix::Zero(1, 1));

        for ( int l = 0; l < s.size(); ++l ) {
            if ( overlapMap[l].size() == 0 ) {
                continue;
            }

            Matrix & o = overlapMap[l][idxBox][c];
            Matrix & loss = lmap[l];

            loss = Matrix::Zero(o.rows(), o.cols());

            // PASCAL VOC loss
            loss = (o.array() < 0.5F).select(1, loss);

            // fg overlap
            // Require at least some overlap with the foreground bounding box
            //    In an image with multiple objects, this constraint encourages a
            //    diverse set of false positives (otherwise, they will tend to come
            //    from the same high-scoring / low-overlapping region of the image
            //    -- i.e. somewhere in the background)
            loss = (o.array() < fgOverlap).select(-Inf, loss);

            // bg overlap
            // Mark root locations that have too much overlap with background boxes as invalid
            //     We don't want to select detections of other foreground objects
            //     in the image as false positives (i.e., no true positive should
            //     be allowed to be used as a false positive)
            for ( int i = 0; i< nbBoxes; ++i ) {
                if ( i == idxBox ) {
                    continue;
                }

                Matrix & obg = overlapMap[l][i][c];

                loss = (obg.array() >= bgOverlap).select(-Inf, loss);
            }

            // loss adjustment
            s[l] += loss;

            if ( needCopy ) {
                setScoreMapStatus(objComp, l);
            }
        } // for l
    } // for c

    computeORNode(grammar_->getRootNode());
}

void Inference::release()
{
    scoreMaps_.clear();
    scoreMapCopies_.clear();
    scoreMapStatus_.clear();
    deformationX_.clear();
    deformationY_.clear();
    lossMaps_.clear();
}


} // namespace RGM
