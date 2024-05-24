#include "header.h"


void initDataSet(DataSet *ds)
{
    memset(ds->n_points, 0, sizeof(ds->n_points));
}

void importData(DataSet *ds, Point *data_array, int data_array_size)
{
    //need data_array_size. we can't eliminate the parameter and count the total number of points because we are adding points, we don't know how many

    for (index_type i = 0; i < data_array_size; i++)
    {
        index_type layer = data_array[i].layer_num - 1;
        ds->array[layer][ds->n_points[layer]+1] = data_array[i];
        ds->n_points[layer]++; //here n_points is not counting the blank spot at index 0. 
        
    }
    // iterating over the layers in DataSet
    for (index_type i = 0; i < num_layers; i++)
    {
        //sorts the points in the ith layer
        qsort(&ds->array[i][1], ds->n_points[i], sizeof(Point), comparePoints);
    }
}

void addBoundaryPoint(DataSet *ds, float offset)
{
    ds->boundaryPoint_offset = offset;

    for (index_type i = 0; i < num_layers; i++) {
        //adding two boundary points in each layer
        // inserting at the beginning
        ds->array[i][0].layer_num = i + 1;
        ds->array[i][0].radius = (i + 1) * 5;
        ds->array[i][0].phi = ds->array[i][1].phi; // 1 gets what was originally the index 0, but is not index 1 due to the insertion.
        ds->array[i][0].z = -1 * ((trapezoid_edges[i]) - offset) - offset; //trapezoid edges is constant and initialized with the offset added. to preserve the original statement, we do it like this

        // appending at the end
        index_type lastIndex = ds->n_points[i] + 1; // after shifting, there's one more point
        ds->array[i][lastIndex].layer_num = i + 1;
        ds->array[i][lastIndex].radius = (i + 1) * 5;
        ds->array[i][lastIndex].phi = ds->array[i][1].phi; // 1 gets what was originally the index 0, but is not index 1 due to the insertion.
        ds->array[i][lastIndex].z = trapezoid_edges[i]; //here we want x.0001

        //now factors in the addition of both boundary points because n_points previously was counting true point additions, and did not count the blank index 0.
        ds->n_points[i] += 2;    

        // adjusting positions using insertion sort techniques as opposed to sorting the entire array. 
        // we have the guarentee from importData that the array was sorted
        // assigned points to indices first to avoid risk of comparing uninitialized "blank" points.
        // as opposed to full sorting algorithms like mergesort, each call here is O(N) and has the potential to escape much earlier. 
        adjustPointPositionFront(ds->array[i], ds->n_points[i], 0); // adjust the start boundary
        adjustPointPositionBack(ds->array[i], ds->n_points[i], lastIndex); // adjust the end boundary
    }

}

void adjustPointPositionFront(Point *array, int n_points, int start_index) {
    // move the point at start_index to its correct position to maintain sorted order
    Point toInsert = array[start_index];
    int j = start_index;
    //by checking if j < n_points-2, we are not going to the last index, thus, we will never have a situation where the end boundary point gets shifted before we position it
    //we check n_points-2 instead of n_points-1 because we have j+1 logic, j is the baseline and the comparison is with the next index, so we need j+1 to be not the end, but 1 away from it.
    //this is a valid cutoff because the z is the primary (first [and only in the case of the implementation, non-debugging comparator]) comparison in the comparator, and the trapezoid edges are always positive integers, so -x < x when x is a positive integer. 
    //it cannot be 0, so there is no possible equality as well, which could affect the debugging comparator.
    while (j < n_points - 2 && comparePoints(&array[j + 1], &toInsert) < 0) { // once we find one element does not need to be moved, we can stop, because the array is monotonic because it is sorted
        array[j] = array[j + 1]; // shift elements left, the other element(s) should come before the boundary point
        j++;
    }
    array[j] = toInsert; // place the element at its correct position
}

void adjustPointPositionBack(Point *array, int n_points, int start_index) {
    // move the point at start_index to its correct position to maintain sorted order
    Point toInsert = array[start_index];
    int j = start_index;
    //similarly, j > 1 ensures it doesn't reach the first index [j will end at 1 after checking if 2 should swap with 1], which while it wouldn't throw off the front position such that the adjustFront method doesn't work because that has already been called,
    //it is beneficial not to check the first index because it is a pointless computation. we can guarentee it will not shift.
    while (j > 1 && comparePoints(&array[j - 1], &toInsert) > 0) { // once we find one element does not need to be moved, we can stop, because the array is monotonic because it is sorted
        array[j] = array[j - 1]; // shift elements right, the other element(s) should come after the boundary point
        j--;
    } 
    array[j] = toInsert; // place the element at its correct position
}


