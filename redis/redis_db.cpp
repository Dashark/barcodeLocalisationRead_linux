#include <assert.h>
#include <hiredis.h>

#include "redis_db.h"


struct redisHost {
    std::string ip;
    int port;
    int db_idx;
};

static redisContext *g_ctx = NULL;
static redisHost *g_addr = NULL;

bool connect(std::string ip, int port, int db_idx) {
    g_addr = new redisHost{ip, port, db_idx};

    struct timeval timeout = { 1, 500000 };
    g_ctx = redisConnectWithTimeout(ip.c_str(), port, timeout);
    if (g_ctx == NULL || g_ctx->err) {
        if (g_ctx) {
            printf("Connection error: %s\n", g_ctx->errstr);
            redisFree(g_ctx);
            g_ctx = NULL;
        } else {
            printf("Connection error: can't allocate redis context\n");
        }
        assert(false);
        return false;
    }
    
    if (db_idx != 0) {
        redisReply *reply;
        reply = (redisReply*) redisCommand(g_ctx,"SELECT %d", db_idx);
        if (NULL == reply) {
            if (g_ctx->err) {
                printf("SELECT error: %s\n", g_ctx->errstr);
            }
            return false;
        } else {
            if (reply->type == REDIS_REPLY_ERROR) {
                printf("SELECT reply error: %s\n", reply->str);
                assert(false);
                return false;
            }
        }
        freeReplyObject(reply);
    }

    return true;
}

void disconnect() {
    if (g_ctx != NULL) {
        redisFree (g_ctx); //Disconnects and frees the context.
        g_ctx = NULL;
    }
}

bool reconnect() {
    disconnect();
    if (NULL != g_addr) {
        return connect(g_addr->ip.c_str(), g_addr->port, g_addr->db_idx);
    } else {
        return false;
    }
}

bool set_bin(std::string key, std::string& value) {
    assert(NULL != g_ctx);
    redisReply *reply;
    reply = (redisReply*) redisCommand(g_ctx,"SET %s %b", key.c_str(), value.data(), value.size());
    if (NULL == reply) {
        if (g_ctx->err) {
            printf("SET error: %s\n", g_ctx->errstr);
        }
        return false;
    } else {
        if (reply->type == REDIS_REPLY_ERROR) {
            printf("SET reply error: %s\n", reply->str);
            assert(false);
            return false;
        }
    }
    freeReplyObject(reply);

    return true;
}

bool setex_bin(std::string key, int seconds, std::string& value) {
    assert(NULL != g_ctx);
    redisReply *reply;
    reply = (redisReply*) redisCommand(g_ctx,"SETEX %s %d %b", key.c_str(), seconds, value.data(), value.size());
    if (NULL == reply) {
        if (g_ctx->err) {
            printf("SETEX error: %s\n", g_ctx->errstr);
        }
        return false;
    } else {
        if (reply->type == REDIS_REPLY_ERROR) {
            printf("SETEX reply error: %s\n", reply->str);
            assert(false);
            return false;
        }
    }
    freeReplyObject(reply);

    return true;
}

std::string get_bin(std::string key) {
    assert(NULL != g_ctx);
    std::string value;
    redisReply *reply;
    reply = (redisReply*) redisCommand(g_ctx,"GET %s", key.c_str());
    if (NULL == reply) {
        if (g_ctx->err) {
            printf("GET error: %s\n", g_ctx->errstr);
        }
        return value;
    } else {
        if (reply->type == REDIS_REPLY_ERROR) {
            printf("GET reply error: %s\n", reply->str);
            assert(false);
            return value;
        } else if (reply->type == REDIS_REPLY_STRING) {
            value = std::string(reply->str, reply->len);
        } else if (reply->type == REDIS_REPLY_NIL) {
            printf("GET reply is NIL\n");
        }
    }
    freeReplyObject(reply);

    return value;
}

bool lpush_bin(std::string key, std::string& value) {
    assert(NULL != g_ctx);
    redisReply *reply;
    reply = (redisReply*) redisCommand(g_ctx,"LPUSH %s %b", key.c_str(), value.data(), value.size());
    if (NULL == reply) {
        if (g_ctx->err) {
            printf("LPUSH error: %s\n", g_ctx->errstr);
        }
        return false;
    } else {
        if (reply->type == REDIS_REPLY_ERROR) {
            printf("LPUSH reply error: %s\n", reply->str);
            assert(false);
            return false;
        }
    }
    freeReplyObject(reply);

    return true;
}

std::string rpop_bin(std::string key) {
    assert(NULL != g_ctx);
    std::string value;
    redisReply *reply;
    reply = (redisReply*) redisCommand(g_ctx,"RPOP %s", key.c_str());
    if (NULL == reply) {
        if (g_ctx->err) {
            printf("RPOP error: %s\n", g_ctx->errstr);
        }
        return value;
    } else {
        if (reply->type == REDIS_REPLY_ERROR) {
            printf("RPOP reply error: %s\n", reply->str);
            assert(false);
            return value;
        } else if (reply->type == REDIS_REPLY_STRING) {
            value = std::string(reply->str, reply->len);
        } else if (reply->type == REDIS_REPLY_NIL) {
            printf("RPOP reply is NIL\n");
        }
    }
    freeReplyObject(reply);

    return value;
}
