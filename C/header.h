// gdb bin/ProcessInput
// run < CPP/wedgeData_v3_128.txt
// bt
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>


//make conversion ratio to micron macro

#ifdef MAIN_C
#define EXTERN
#else
#define EXTERN extern
#endif

#define min(X, Y) ((X) < (Y) ? (X) : (Y))
#define max(X, Y) ((X) < (Y) ? (Y) : (X))

#define index_type int //change to unsigned int once it is verified that there are no errors
#define CONVERSION_TYPE long

#define MAX_LAYERS 5
#define MAX_POINTS_IN_EVENT 512
#define MAX_POINTS_PER_LAYER 256    // max size of vector of points "vect" in CPP. equivalent to MAX_POINTS_PER_DATASET
#define MAX_POINTS_FOR_DATASET MAX_POINTS_PER_LAYER    // max size of vector of points "vect" in CPP

#define MAX_POINTS_IN_LINE MAX_LAYERS // a point on the line is calculated for each layer in the environment.
#define MAX_POINTS_IN_SUPERPOINT 16
#define MAX_SUPERPOINTS_IN_PATCH 5
#define MAX_PARALLELOGRAMS_PER_PATCH MAX_LAYERS - 1 // layer 1 is a vertical ribbon, the other 4 layers are sloping, so each intersects with layer 1 to make a parallelogram
#define MAX_PATCHES 32                              // upper bound, 14-18 average.
// #define MAX_LINES __ //only used in visualization
#define MAX_SUPERPOINTS_IN_COVER (MAX_PATCHES * MAX_SUPERPOINTS_IN_PATCH)

#ifdef MAIN_C
    const index_type CONVERSION_FACTOR = 1000000; //10k = 39k diff, 100k = 45k diff
    const CONVERSION_TYPE top_layer_lim = 50 * CONVERSION_FACTOR;
    const CONVERSION_TYPE beam_axis_lim = 15 * CONVERSION_FACTOR;
    const index_type num_layers = 5;
    const index_type radii[MAX_LAYERS] = {
        5 * CONVERSION_FACTOR, 
        10 * CONVERSION_FACTOR, 
        15 * CONVERSION_FACTOR, 
        20 * CONVERSION_FACTOR, 
        25 * CONVERSION_FACTOR
    };
    const CONVERSION_TYPE parallelogramSlopes[MAX_LAYERS-1] = {
        0,
        (radii[0] - radii[1]) / (radii[4] - radii[1]), 
        (radii[0] - radii[2]) / (radii[4] - radii[2]), 
        (radii[0] - radii[3]) / (radii[4] - radii[3])
    };
    const CONVERSION_TYPE radii_leverArm[MAX_LAYERS-1] = {
        1 - parallelogramSlopes[0], 
        1 - parallelogramSlopes[1], 
        1 - parallelogramSlopes[2], 
        1 - parallelogramSlopes[3]
    };
    const CONVERSION_TYPE trapezoid_edges[MAX_LAYERS] = {
        22.0001 * CONVERSION_FACTOR,
        29.0001 * CONVERSION_FACTOR,
        36.0001 * CONVERSION_FACTOR,
        43.0001 * CONVERSION_FACTOR,
        50.0001 * CONVERSION_FACTOR
        //radii[0] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR, 
        //radii[1] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR, 
        //radii[2] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR, 
        //radii[3] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR, 
        //radii[4] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR
    };
#else
    extern index_type radii[MAX_LAYERS];
    extern CONVERSION_TYPE parallelogramSlopes[MAX_LAYERS-1];
    extern CONVERSION_TYPE radii_leverArm[MAX_LAYERS-1];
    extern CONVERSION_TYPE trapezoid_edges[MAX_LAYERS];
    extern index_type CONVERSION_FACTOR;
    extern CONVERSION_TYPE top_layer_lim;
    extern CONVERSION_TYPE beam_axis_lim;
    extern index_type num_layers;
#endif


typedef struct
{
    index_type layer_num;
    index_type radius;
    CONVERSION_TYPE phi;
    CONVERSION_TYPE z;
} Point;

typedef struct
{
    Point array[MAX_LAYERS][MAX_POINTS_FOR_DATASET]; // 2D array of points
    int n_points[MAX_LAYERS];                        // number of points in each layer of the array
    //index_type total_points; //not used
    CONVERSION_TYPE boundaryPoint_offset;
} DataSet;

typedef struct
{
    index_type layer_num;
    CONVERSION_TYPE pSlope;

    CONVERSION_TYPE shadow_bottomL_jR;
    CONVERSION_TYPE shadow_bottomR_jR;
    CONVERSION_TYPE shadow_bottomL_jL;
    CONVERSION_TYPE shadow_bottomR_jL;

    CONVERSION_TYPE z1_min;
    CONVERSION_TYPE z1_max;
} Parallelogram;

typedef struct
{
    Point points[MAX_POINTS_IN_SUPERPOINT];
    CONVERSION_TYPE z_values[MAX_POINTS_IN_SUPERPOINT];
    index_type point_count;
    CONVERSION_TYPE min;
    CONVERSION_TYPE max;
} wedgeSuperPoint;

typedef struct
{
    int end_layer;
    int left_end_layer;
    int right_end_layer;
    CONVERSION_TYPE left_end_lambdaZ;
    CONVERSION_TYPE right_end_lambdaZ;
    CONVERSION_TYPE apexZ0;

    CONVERSION_TYPE shadow_fromTopToInnermost_topL_jL;
    CONVERSION_TYPE shadow_fromTopToInnermost_topL_jR;
    CONVERSION_TYPE shadow_fromTopToInnermost_topR_jL;
    CONVERSION_TYPE shadow_fromTopToInnermost_topR_jR;

    CONVERSION_TYPE a_corner[2];
    CONVERSION_TYPE b_corner[2];
    CONVERSION_TYPE c_corner[2];
    CONVERSION_TYPE d_corner[2];

    wedgeSuperPoint superpoints[MAX_SUPERPOINTS_IN_PATCH]; //changed to direct assignment as opposed to pointer
    index_type superpoint_count;

    bool flatBottom;
    bool flatTop;

    bool squareAcceptance;
    bool triangleAcceptance;

    Parallelogram parallelograms[MAX_PARALLELOGRAMS_PER_PATCH];
    index_type parallelogram_count;
} wedgePatch;

extern int Point_load(Point *p);
extern void importData();
extern void addBoundaryPoint(CONVERSION_TYPE offset);
extern void initWedgeSuperPoint(wedgeSuperPoint *wsp, Point *points, int pointCount);
extern int areWedgeSuperPointsEqual(wedgeSuperPoint *wsp1, wedgeSuperPoint *wsp2);
extern void wedgePatch_init(wedgePatch *wp, wedgeSuperPoint *superpointsI, int superpoint_count, CONVERSION_TYPE apexZ0I);
extern CONVERSION_TYPE straightLineProjectorFromLayerIJtoK(CONVERSION_TYPE z_i, CONVERSION_TYPE z_j, int i, int j, int k);
extern CONVERSION_TYPE straightLineProjector(CONVERSION_TYPE z_top, CONVERSION_TYPE z_j, int j);
extern void getParallelograms(wedgePatch *wp);
extern void getShadows(wedgePatch *wp, CONVERSION_TYPE zTopMin, CONVERSION_TYPE zTopMax);
extern void get_acceptanceCorners(wedgePatch *wp);
extern void get_end_layer(wedgePatch *wp);
extern void initWedgeCover();
extern int comparePoints(const void *a, const void *b);
extern void add_patch(wedgePatch *curr_patch);
extern void delete_patch(int index);
extern index_type get_index_from_z(int layer, CONVERSION_TYPE z_value);
extern void solve(CONVERSION_TYPE apexZ0, int ppl, bool leftRight);
extern void makePatches_ShadowQuilt_fromEdges(int ppl, bool leftRight);
extern void makePatch_alignedToLine(CONVERSION_TYPE apexZ0, CONVERSION_TYPE z_top, int ppl, bool leftRight, bool CONVERSION_TYPE_middleLayers_ppl);
extern void wedge_test(CONVERSION_TYPE apexZ0, int ppl, int wedges[]);
extern int CONVERSION_TYPECompare(const void *a, const void *b);

EXTERN DataSet Gdata;
EXTERN wedgePatch patches[MAX_PATCHES];
EXTERN index_type n_patches;
