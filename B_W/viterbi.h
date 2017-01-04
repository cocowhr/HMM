#include <WinSock2.h>  
#include "mysql.h"  
#pragma comment(lib,"wsock32.lib")  
#pragma comment(lib,"libmysql.lib")   
#include <iostream>  
#include <vector>  
#include <map>  
#include <string>  
#include <algorithm>
#include <atlstr.h>
using namespace std;
#define ITEMNUM 1285
#define ITEMLEN 4
#define SEQUENCELEN 4
#define TV 0.5	//短序列匹配阈值
typedef pair<int,double> seq;
/***mysql变量***/
MYSQL mysql;  
MYSQL_RES *result;  
MYSQL_ROW row;  
CString strSQL;
int res;