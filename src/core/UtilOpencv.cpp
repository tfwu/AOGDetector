#include <boost/math/special_functions.hpp>

#include "UtilOpencv.hpp"


namespace RGM
{

// -------- OpencvUtil -------

cv::Mat OpencvUtil::pictureHOG(cv::Mat_<Scalar> & filter, int bs)
{
    // construct a "glyph" for each orientaion
    cv::Mat_<Scalar> bim1(bs, bs, Scalar(0));

    int hbs = floor(bs/2.0f + 0.5f);
    cv::Range rg(hbs-1, hbs+1);
    bim1.colRange(rg) = 1;

    int ht     = filter.size[0];
    int wd     = filter.size[1];
    int numOri = filter.size[2];

    float degreeUnit = (float)180.0F / numOri;

    std::vector<cv::Mat_<Scalar> > bim(numOri);
    bim[0] = bim1;

    for (int i=1; i<numOri; ++i) {
        bim[i] = OpencvUtil::rotate(cv::Mat(bim1), -i*degreeUnit);
    }

    // make pictures of positive weights bs adding up weighted glyphs
    cv::Mat_<Scalar> posf = cv::max(filter, 0);
    cv::Mat_<Scalar> img(bs*ht, bs*wd, Scalar(0));
    Scalar maxV = -1;
    for ( int i=0; i<ht; ++i ) {
        for ( int j=0; j<wd; ++j ) {
            for ( int k=0; k<numOri; ++k) {
                img(cv::Rect(j*bs, i*bs, bs, bs)) += (bim[k] * posf(i, j, k));
                maxV = std::max<Scalar>(maxV, posf(i, j, k));
            }
        }
    }

    img *= (255.0F / maxV);

    /*cv::Mat imgShow;
    cv::normalize(img, imgShow, 255, 0.0, CV_MINMAX, CV_8UC1);

    cv::imshow("Debug", imgShow);
    cv::waitKey(0);*/

    return cv::Mat(img);
}

cv::Mat OpencvUtil::rotate(cv::Mat m, float degree)
{
    cv::Point2f src_center(m.cols/2.0F, m.rows/2.0F);
    cv::Mat rot_mat = cv::getRotationMatrix2D(src_center, degree, 1.0);
    cv::Mat dst;
    cv::warpAffine(m, dst, rot_mat, m.size());
    return dst;
}


cv::Mat OpencvUtil::subarray(cv::Mat img, cv::Rect roi, float padFactor, int padType) {

    int ht = ROUND( roi.height * (1+padFactor) );
    int wd = ROUND( roi.width * (1+padFactor) );

    int padx  = (wd - roi.width) / 2;
    int padx1 = wd - roi.width - padx;
    int pady  = (ht - roi.height) / 2;
    int pady1 = ht - roi.height - pady;

    int x1 = roi.x - padx;
    int x2 = roi.br().x + padx1;
    int y1 = roi.y - pady;
    int y2 = roi.br().y + pady1;

    cv::Rect newroi(x1, y1, x2-x1, y2-y1);

    return subarray(img, newroi, padType);
}

cv::Mat OpencvUtil::subarray(cv::Mat img, cv::Rect roi, int padType)
{
    int x1 = roi.x;
    int x2 = roi.br().x;
    int y1 = roi.y;
    int y2 = roi.br().y;

    if (x1>=0 && y1>=0 && x2<=img.cols && y2<=img.rows) {
        return img(cv::Rect(x1, y1, x2-x1, y2-y1)).clone();
    } else {
        cv::Mat subImg(y2-y1, x2-x1, img.type(), cv::Scalar::all(0));

        if (padType==1) {
            for ( int y=y1; y<y2; ++y ) {
                int yy = std::min<int>(img.rows-1, std::max<int>(0, y));
                for ( int x=x1; x<x2; ++x ) {
                    int xx = std::min<int>(img.cols-1, std::max<int>(0, x));

                    if (img.channels()==3) {
                        subImg.at<cv::Vec3b>(y-y1, x-x1) = img.at<cv::Vec3b>(yy, xx);
                    } else {
                        subImg.at<unsigned char>(y-y1, x-x1) = img.at<unsigned char>(yy, xx);
                    }
                }
            }
        } else {
            int xx1 = std::max<int>(0, x1);
            int xx2 = std::min<int>(img.cols, x2);
            int yy1 = std::max<int>(0, y1);
            int yy2 = std::min<int>(img.rows, y2);
            int wd = xx2 - xx1;
            int ht = yy2 - yy1;

            int padx = xx1 - x1;
            int pady = yy1 - y1;

            img(cv::Rect(xx2, yy1, wd, ht)).copyTo( subImg(cv::Rect(padx, pady, wd, ht)) );
        }

        return subImg;
    }
}


// -------- OpencvUtil_ -------

template<typename T>
cv::Mat_<T> OpencvUtil_<T>::resize(cv::Mat_<T> & m, float factor, int method)
{
    int ht = floor(m.size[0] * factor + 0.5f);
    int wd = floor(m.size[1] * factor + 0.5f);

    int dims[] = {ht, wd, m.size[2]};

    cv::Mat_<T> result(3, dims);

    for ( int d=0; d<dims[2]; ++d ) {
        cv::Mat m1(m.size[0], m.size[1], CV_MAKE_TYPE(cv::DataDepth<T>::value, 1));
        for ( int r=0; r<m1.rows; ++r) {
            for ( int c=0; c<m1.cols; ++c ) {
                m1.at<T>(r, c) = m(r, c, d);
            }
        }

        cv::Mat rm;
        cv::resize(m1, rm, cv::Size(wd, ht), 0, 0, method);

        for ( int r=0; r<ht; ++r) {
            for ( int c=0; c<wd; ++c) {
                result(r, c, d) = rm.at<T>(r, c);
            }
        }
    }

    return result;
}

// instantiation
template class OpencvUtil_<int>;
template class OpencvUtil_<float>;
template class OpencvUtil_<double>;

} // namespace RGM
