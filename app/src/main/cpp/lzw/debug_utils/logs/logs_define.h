#if defined(MY_DEBUG)
#define LOG_DEBUG(args...) \
        write_logs (args); \

#else
#define LOG_DEBUG(args...);
#endif

//-----------------------------------------------------------------------------------------
