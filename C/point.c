#include "header.h"

int Point_load(Point* p)
{
	//reads input of the form (layer_num,radius,phi,z) to populate point structure. 1 if worked, 0 if not.
	if (scanf("(%d,%f,%f,%f)", &p->layer_num, &p->radius, &p->phi, &p->z) == 4) //p->layer_num equivalent to (*p).layer_num, address of that
	{
		return 1;
	}

	return 0;
}

int comparePoints(const void* a, const void* b) {
    const Point* pointA = (const Point*)a;
    const Point* pointB = (const Point*)b;
    if (pointA->z < pointB->z) return -1;
    if (pointA->z > pointB->z) return 1;
    return 0;
}
