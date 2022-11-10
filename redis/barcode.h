#ifndef _BARCODE_H_
#define _BARCODE_H_

#include <opencv2/opencv.hpp>
#include <string>

typedef struct {
  std::string lhs, rhs;
  cv::Rect crop;
} BAR_CODE;

BAR_CODE barcode(cv::Mat &img);

#endif /* _BARCODE_H_ */