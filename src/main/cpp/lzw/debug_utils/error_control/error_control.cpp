#include "error_control.hpp"

//-----------------------------------------------------------------------------------------

int is_nullptr (const void* ptr) {
    if (ptr == nullptr) {
        return 1;
    }
    else return 0;
}

int print_error_message (const char* file_, const char* func_, int line_) {
    const char* error_mes = "\033[91mERROR Message\033[0m";
    fprintf (stderr, "-------------------------%s----------------------\n", error_mes);
    fprintf (stderr, "| filename:     %s %s %s\n", BLUE_C,   file_, RESET_C);
    fprintf (stderr, "| name_of_func: %s %s %s\n", GREEN_C,  func_, RESET_C);
    fprintf (stderr, "| num_of_line:  %s %d %s\n", YELLOW_C, line_, RESET_C);
    fprintf (stderr, "------------------------------------------------------------\n\n");

    return 0;
}

//-----------------------------------------------------------------------------------------
