#include "messager.h"
#include "redis_db.h"

static bool redis_connected = false;
static std::string redis_ip = "127.0.0.1";
static int redis_port = 6379, redis_db = 1;

static std::string request_key = "barcode:request";
static std::string response_key = "barcode:response";

bool fetch_request(cv::Mat& image) {
    if (!redis_connected) {
        redis_connected = connect(redis_ip, redis_port, redis_db);
    }

    if (redis_connected) {
        std::string img_key = rpop_bin(request_key);
        if (img_key.empty()) {
            return false;
        }
        std::string str_decode = get_bin(img_key);
        if (str_decode.empty()) {
            return false;
        }
        std::vector<uchar> data_decode(str_decode.begin(), str_decode.end());
        cv::Mat img_decode = cv::imdecode(data_decode, cv::IMREAD_COLOR);
        image = img_decode;
        return true;
    } else {
        return false;
    }
}

bool push_response(std::string message) {
    if (!redis_connected) {
        redis_connected = connect(redis_ip, redis_port, redis_db);
    }

    if (redis_connected) {
        if (lpush_bin(response_key, message)) {
            return true;
        } else {
            redis_connected = reconnect();
            return false;
        }
    } else {
        return false;
    }
}

