#ifndef ERROR_CONTR_H
#define ERROR_CONTR_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <chrono>

#include "console_colors.hpp"
#include "ASSERT.hpp"

//-----------------------------------------------------------------------------------------

#define CUR_POS_IN_PROG __FILE__, __PRETTY_FUNCTION__, __LINE__

//-----------------------------------------------------------------------------------------

int   is_nullptr (const void* ptr);
int   print_error_message (const char* file_, const char* func_, int line_);

#endif
