#include "mongoose.h"
#include <opencv2/opencv.hpp>

#include "../redis/barcode.h"

static const char *s_http_addr = "http://0.0.0.0:8080";    // HTTP port


static void fn(struct mg_connection *c, int ev, void *ev_data, void *fn_data) {
    if (MG_EV_POLL == ev) {
        return;
    }
  
    if (ev == MG_EV_HTTP_MSG) {
        struct mg_http_message *hm = (struct mg_http_message *) ev_data;
#if 0
        printf("method: %.*s\n", (int) hm->method.len, hm->method.ptr);
        printf("uri: %.*s\n", (int) hm->uri.len, hm->uri.ptr);
        printf("body: %.*s\n", (int) hm->body.len, hm->body.ptr);
#endif
        if (mg_http_match_uri(hm, "/api")) {
            long height = mg_json_get_long(hm->body, "$.height", -1);
            long width = mg_json_get_long(hm->body, "$.width", -1);
            
            int len = -1;
            char* data = mg_json_get_b64(hm->body, "$.data", &len); // base64 decoded
            assert(len == height * width * 3);
            printf("height: %ld, width: %ld, len: %d\n", height, width, len);

            if (data) {
                printf("len: %d\n", len);
                cv::Mat frame(height, width, CV_8UC3, data);
                BAR_CODE bc = barcode(frame);
                printf("%s %s", bc.lhs, bc.rhs);
                cv::imshow("cpp", frame);
                cv::waitKey();
                cv::destroyAllWindows();
            } else {
                printf("data is NUL\n");
            }

            if (data) {
                free(data);
            }
            
            mg_http_reply(c, 200, "Content-Type: application/json\r\n", 
                            "{\"result\": \"%.*s\"}\n", (int) hm->uri.len, hm->uri.ptr);
        } else {
            mg_http_reply(c, 404, "", "%s", "Not Found\n");
        }
    }

    (void) fn_data;
}

int main(void) {
    struct mg_mgr mgr;                            // Event manager
    mg_log_set(MG_LL_INFO);                       // Set log level
    mg_mgr_init(&mgr);                            // Initialise event manager
    mg_http_listen(&mgr, s_http_addr, fn, NULL);  // Create HTTP listener
    for (;;) mg_mgr_poll(&mgr, 1000);             // Infinite event loop
    mg_mgr_free(&mgr);
    return 0;
}
