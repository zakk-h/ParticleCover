#include "header.h"

void wedge_test(CONVERSION_TYPE apexZ0, int ppl, int wedges[])
{
    int numEventsLoaded = 0;

    FILE *myfile;
    myfile = fopen("C/cOutput.txt", "w"); 

    if (myfile == NULL)
    {
        printf("Error opening file");
        return;
    }

    for (index_type z = 0; z < wedges[1]; ++z)
    { 
        printf("wedge %d\n", z); //main print
        fprintf(myfile, "wedge %d\n", z); //file to diff

        importData();

        if(z<wedges[0]) continue;
        
        addBoundaryPoint(0.0001*CONVERSION_FACTOR); // with default param

        initWedgeCover();

        solve(apexZ0, ppl, false); // solve modifies  false is from the left right align (previously a parameter in wedge test)

        for (int i = 0; i < n_patches; ++i)
        {
            fprintf(myfile, "Patch \n");
            fprintf(myfile, "%ld\n", lround((patches[i].shadow_fromTopToInnermost_topL_jL * 10000)/CONVERSION_FACTOR));
            fprintf(myfile, "%ld\n", lround((patches[i].shadow_fromTopToInnermost_topL_jR * 10000)/CONVERSION_FACTOR));
            fprintf(myfile, "%ld\n", lround((patches[i].shadow_fromTopToInnermost_topR_jL * 10000)/CONVERSION_FACTOR));
            fprintf(myfile, "%ld\n", lround((patches[i].shadow_fromTopToInnermost_topR_jR * 10000)/CONVERSION_FACTOR));

            for (int j = 0; j < patches[i].superpoint_count; ++j)
            {
                fprintf(myfile, "Superpoint \n");
                for (int r = 0; r < patches[i].superpoints[j].point_count; ++r)
                {
                    Point currentPt = patches[i].superpoints[j].points[r];
                    fprintf(myfile, "%d %.4f %d %.4f\n",
                            currentPt.layer_num,
                            (float)currentPt.phi/CONVERSION_FACTOR,
                            (float)currentPt.radius/CONVERSION_FACTOR,
                            (float)currentPt.z/CONVERSION_FACTOR);
                }
            }
        }
        for (int i = 0; i < n_patches; ++i)
        {
            fprintf(myfile, "[%ld, %ld]\n",
                    lround(patches[i].a_corner[0] * 10000/CONVERSION_FACTOR),
                    lround(patches[i].a_corner[1] * 10000/CONVERSION_FACTOR));
            fprintf(myfile, "[%ld, %ld]\n",
                    lround(patches[i].b_corner[0] * 10000/CONVERSION_FACTOR),
                    lround(patches[i].b_corner[1] * 10000/CONVERSION_FACTOR));
            fprintf(myfile, "[%ld, %ld]\n",
                    lround(patches[i].c_corner[0] * 10000/CONVERSION_FACTOR),
                    lround(patches[i].c_corner[1] * 10000/CONVERSION_FACTOR));
            fprintf(myfile, "[%ld, %ld]\n",
                    lround(patches[i].d_corner[0] * 10000/CONVERSION_FACTOR),
                    lround(patches[i].d_corner[1] * 10000/CONVERSION_FACTOR));
            fprintf(myfile, "\n");
        }
        // instead of making an array of all events and passing them in, we only need access to them individually, so we will loop through and process as we create them.
    }

    fclose(myfile);
}
