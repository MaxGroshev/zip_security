#pragma once

//-----------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>

#include "logs_define.h"

//-----------------------------------------------------------------------------------------

template <class... Args>
int   write_logs (Args... log_text);

#include "logs.tpp"
