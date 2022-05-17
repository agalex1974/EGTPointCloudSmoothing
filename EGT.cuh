//
// Created by agalex on 11/10/21.
//

#ifndef SMOOTHINGOFPOINTCLOUD_EGT_CUH
#define SMOOTHINGOFPOINTCLOUD_EGT_CUH

void EGTsmoothing(float* in_x, float* in_y, float* in_z, int count, float lambda, float mu,
                  int* neighbors, int maxNeighbors, float* out_x, float* out_y, float* out_z,
                  int iterationCount, int isRegularized, float ratio);

#endif //SMOOTHINGOFPOINTCLOUD_EGT_CUH
