#include "header.h"

void initWedgeSuperPoint(wedgeSuperPoint *wsp, Point *points, int pointCount)
{
    wsp->point_count = pointCount;
    wsp->min = LONG_MAX; //assigned opposite values so they get replaced instantly
    wsp->max = -LONG_MAX;

    // more efficient approach
    // instead of making a temp array and then transferring contents, add the values directly
    // simultaneously, you can determine the min and max, as opposed to doing it after the fact
    for (int i = 0; i < pointCount; i++)
    {
        wsp->points[i] = points[i];
        wsp->z_values[i] = points[i].z;

        if (points[i].z < wsp->min)
            wsp->min = points[i].z;
        if (points[i].z > wsp->max)
            wsp->max = points[i].z;
    }
}

// operator overloading not allowed in C, write separate method to check equality
int areWedgeSuperPointsEqual(wedgeSuperPoint *wsp1, wedgeSuperPoint *wsp2)
{
    //return (wsp1->min == wsp2->min) && (wsp1->max == wsp2->max);
    const CONVERSION_TYPE tolerance = 0.0001*CONVERSION_FACTOR;
    return (fabs(wsp1->min - wsp2->min) < tolerance) && (fabs(wsp1->max - wsp2->max) < tolerance);
}
