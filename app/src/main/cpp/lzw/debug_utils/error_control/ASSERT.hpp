#pragma once
#if defined(MY_DEBUG)
#define ASSERT(condition)                         \
;{                                                \
        if (!(condition))                         \
        {                                         \
            print_error_message (CUR_POS_IN_PROG);\
            exit (1);                             \
        }                                         \
}                                                 \

#else
#define ASSERT(condition);

#endif
