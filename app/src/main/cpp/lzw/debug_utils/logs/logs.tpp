#ifndef LOGS_TPP
#define LOGS_TPP

//-----------------------------------------------------------------------------------------

template <class... Args>
int write_logs (Args... log_text) {
    //std::cout << std::vformat(std::make_format_args(log_text...));
    (std::cout << ... << log_text);
    return 0;
}

#endif
