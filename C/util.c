#include "header.h"

int CONVERSION_TYPECompare(const void *a, const void *b)
{
    CONVERSION_TYPE diff = *(const CONVERSION_TYPE *)a - *(const CONVERSION_TYPE *)b;
    if (diff < 0)
    {
        return -1;
    }
    else if (diff > 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
