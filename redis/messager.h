#ifndef _MESSAGER_H_
#define _MESSAGER_H_

#include <string>
#include <opencv2/opencv.hpp>

bool fetch_request(cv::Mat& image);
bool push_response(std::string message);

#endif /* _MESSAGER_H_ */
