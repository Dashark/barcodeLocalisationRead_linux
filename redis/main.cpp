#include <unistd.h>
#include "messager.h"
#include "redis_db.h"
#include "barcode.h"

typedef unsigned char byte;

int main(int argc, char* argv[]) {
    for (;;) {
        cv::Mat image;
        if (fetch_request(image)) {
            barcode(image);
            cv::imshow("cpp", image);
            cv::waitKey();
            cv::destroyAllWindows();
        } else {
            usleep(1000 * 1000);
        }
    }

    return 0;
}
