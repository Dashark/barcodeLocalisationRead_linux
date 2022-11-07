#ifndef _HI_REDIS_DB_H_
#define _HI_REDIS_DB_H_

#include <string>

bool connect(std::string ip, int port, int db_idx);
void disconnect();
bool reconnect();

bool set_bin(std::string key, std::string& value);
bool setex_bin(std::string key, int seconds, std::string& value);

std::string get_bin(std::string key);

bool lpush_bin(std::string key, std::string& value);
std::string rpop_bin(std::string key);

#endif /*_HI_REDIS_DB_H_*/