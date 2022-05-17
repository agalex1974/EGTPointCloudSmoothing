#include "CCreoPointCloud.h"
#include <stdlib.h>

int main(int argc, char* argv[])
{
    // Number of neighbors
	int k = atoi(argv[2]);
    // Lambda parameter
	float lambda = atof(argv[3]);
    // The mu parameter 
	float mu = atof(argv[4]);
    // The alpha parameter
	float alpha = atof(argv[5]);
	// 0 for Gaussian weights, 1 if you want point regularization but more smoothing.
    // In the paper it was set to 0.
	int isRegularized = atoi(argv[6]);
    string fileName = "smoothed.txt";
	if (argc > 6) fileName = string(argv[7]);
    printf("Reading point cloud...");
    CCreoPointCloud pc(argv[1], false);
    printf("done.\n");
    printf("smoothing...");
    pc.smooth_iterative(k, lambda, mu, alpha, isRegularized);
    printf("done.\n");
    pc.savePointCloud(fileName.c_str(), false);
    return 0;
}