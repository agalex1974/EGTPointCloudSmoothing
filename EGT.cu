#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <algorithm>

__device__
void normalize(float& x, float& y, float& z){
    float norm = sqrtf(x * x + y * y + z * z);
    x /= norm; y /= norm; z /= norm;
}

__device__
void GetNormalizedPerpendicularVectorToVector(const float& x1, const float& y1, const float& z1,
                                              float& x2, float& y2, float& z2){
    float max = fabs(x1);

    int cordIndex = 0;
    if (max < fabs(y1))
    {
        cordIndex = 1;
        max = fabs(y1);
    }

    if (max < fabs(z1))
    {
        cordIndex = 2;
    }

    x2 = 1.0;
    y2 = 1.0;
    z2 = 1.0;

    switch (cordIndex)
    {
        case 0:
            x2 = (-y1 * y2 - z1 * z2) / x1;
            break;
        case 1:
            y2 = (-x1 * x2 - z1 * z2) / y1;
            break;
        case 2:
            z2 = (-x1 * x2 - y1 * y2) / z1;
            break;
    }
    normalize(x2, y2, z2);
}

__device__
float norm2(const float& x1, const float& y1, const float& z1, const float& x2, const float& y2, const float& z2){
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
}
__device__
void cross(const float& u1, const float& u2, const float& u3,
           const float& v1, const float& v2, const float& v3,
           float& x, float&y, float& z){
    x = u2 * v3 - v2 * u3;
    y = v1 * u3 - u1 * v3;
    z = u1 * v2 - v1 * u2;
}

__device__
void matrix_multiplication(float& a11, float& a12, float& a13, float& a14,
                           float& a21, float& a22, float& a23, float& a24,
                           float& a31, float& a32, float& a33, float& a34,
                           float& x, float& y, float& z){
    float x1 = x; float y1 = y; float z1 = z;
    x = a11 * x1 + a12 * y1 + a13 * z1 + a14;
    y = a21 * x1 + a22 * y1 + a23 * z1 + a24;
    z = a31 * x1 + a32 * y1 + a33 * z1 + a34;
}

__device__
void CreateLocalCoordinateSystem(float& xo, float& yo, float& zo,
                                 const float& xd, const float& yd, const float& zd,
                                 float& a11, float& a12, float& a13, float& a14,
                                 float& a21, float& a22, float& a23, float& a24,
                                 float& a31, float& a32, float& a33, float& a34)
{
    GetNormalizedPerpendicularVectorToVector(xd, yd, zd, a21, a22, a23);
    cross(xd, yd, zd, a21, a22, a23, a31, a32, a33);
    normalize(a31, a32, a33);
    a14 = 0.0;
    a24 = 0.0;
    a34 = 0.0;
    matrix_multiplication(a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, xo, yo, zo);
    a14 = -xo;
    a24 = -yo;
    a34 = -zo;
}

__device__
bool isEllipicGabrielNeighbor(int pnt_idx, int i, float* pnts_x, float* pnts_y, float* pnts_z,
                              const float& x, const float& y, const float& z, int* NNs,
                              int count, float a)
{
    int neigh = NNs[i * count + pnt_idx];
    float xi = pnts_x[neigh]; float yi = pnts_y[neigh]; float zi = pnts_z[neigh];
    float xo = 0.5f * (xi + x); float yo = 0.5f * (yi + y); float zo = 0.5f * (zi + z);

    float d = sqrtf(norm2(xi, yi, zi, x, y, z)) / 2.0f;
    float xaxis_x = xi - x; float xaxis_y = yi - y; float xaxis_z = zi - z;
    normalize(xaxis_x, xaxis_y, xaxis_z);
    float a11, a12, a13, a14;
    float a21, a22, a23, a24;
    float a31, a32, a33, a34;
    a11 = xaxis_x; a12 = xaxis_y; a13 = xaxis_z;
    CreateLocalCoordinateSystem(xo, yo, zo, xaxis_x, xaxis_y, xaxis_z, a11, a12, a13, a14,
                                a21, a22, a23, a24, a31, a32, a33, a34);
    for (int j = 0; j < i; j++)
    {
        neigh = NNs[j * count + pnt_idx];
        xi = pnts_x[neigh]; yi = pnts_y[neigh]; zi = pnts_z[neigh];
        matrix_multiplication(a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, xi, yi, zi);
        float ellipsoidValue = xi * xi + yi * yi / (a * a) + zi * zi / (a * a);
        if (ellipsoidValue < d * d) return false;
    }
    return true;
}

__global__
void calculateEGG(float* pnts_x, float* pnts_y, float* pnts_z, int* NNs, float ratio, int count, int neighborsCount, int start, int batchCount)
{
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t < batchCount){
        int i = t + start;
        float x = pnts_x[i], y = pnts_y[i], z = pnts_z[i];
        for (int j = neighborsCount - 1; j >= 0; j--)
        {
            if (!isEllipicGabrielNeighbor(i, j, pnts_x, pnts_y, pnts_z, x, y, z, NNs, count, ratio))
                NNs[j * count + i] = -1;
        }
    }
}

__global__
void taubin_step(float* in_x, float* in_y, float* in_z,
                 float* out_x, float* out_y, float* out_z, int count, float scale, int* neighbors, int max_neighbors, int isRegularized, int start, int batchCount){
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t < batchCount){
        int i = t + start;
        float cog_x = 0.0f;
        float cog_y = 0.0f;
        float cog_z = 0.0f;
        float sum = 0.0f;
        float x1 = in_x[i]; float y1 = in_y[i]; float z1 = in_z[i];
        for (int n = 0; n < max_neighbors; n++){
            int neigh = neighbors[n * count + i];
            if (neigh != -1) {
                float x2 = in_x[neigh];
                float y2 = in_y[neigh];
                float z2 = in_z[neigh];
                float distance = norm2(x1, y1, z1, x2, y2, z2);
                float w;
                if (isRegularized) w = scale < 0.0f ? 1.0f / (distance + 1e-8f) : 1.0f;
                else w = expf(-distance);
                cog_x += w * (x2 - x1);
                cog_y += w * (y2 - y1);
                cog_z += w * (z2 - z1);
                sum += w;
            }
        }
        out_x[i] = x1 + scale * cog_x / sum;
        out_y[i] = y1 + scale * cog_y / sum;
        out_z[i] = z1 + scale * cog_z / sum;
    }
}

__global__
void scale_points_to_unity(float* x, float* y, float* z, float min, float max, int pointCount)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i < pointCount) {
        float d = max - min;
        x[i] = (x[i] - min) / d;
        y[i] = (y[i] - min) / d;
        z[i] = (z[i] - min) / d;
    }
}

__global__
void produce_output(float* x, float* y, float* z, float* xout, float* yout, float* zout, float min, float max, int pointCount)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i < pointCount) {
        float d = max - min;
        xout[i] = x[i] * d + min;
        yout[i] = y[i] * d + min;
        zout[i] = z[i] * d + min;
    }
}

void EGTsmoothing(float* in_x, float* in_y, float* in_z, int count, float lambda, float mu,
                  int* neighbors, int maxNeighbors, float* out_x, float* out_y, float* out_z,
                  int iterationCount, int isRegularized, float ratio) {
    size_t size = count * sizeof(float);
	size_t sizeNeighbors = count * maxNeighbors * sizeof(int);
    float *din_x, *din_y, *din_z;
    float *dout_x, *dout_y, *dout_z;
    int *dneighbors;
    cudaMalloc((void **) &din_x, size);
    cudaMemcpy(din_x, in_x, size, cudaMemcpyHostToDevice);
    cudaMalloc((void **) &dout_x, size);
    cudaMalloc((void **) &din_y, size);
    cudaMemcpy(din_y, in_y, size, cudaMemcpyHostToDevice);
    cudaMalloc((void **) &dout_y, size);
    cudaMalloc((void **) &din_z, size);
    cudaMemcpy(din_z, in_z, size, cudaMemcpyHostToDevice);
    cudaMalloc((void **) &dout_z, size);
    cudaMalloc((void **) &dneighbors, sizeNeighbors);
    cudaMemcpy(dneighbors, neighbors, sizeNeighbors, cudaMemcpyHostToDevice);
    int max_threads = 2000000;

    float xmin = thrust::reduce(thrust::device, din_x, din_x + count, 1e8f, thrust::minimum<float>());
    float xmax = thrust::reduce(thrust::device, din_x, din_x + count, -1e8f, thrust::maximum<float>());
    float min = xmin;
    float max = xmax;
    float ymin = thrust::reduce(thrust::device, din_y, din_y + count, 1e8f, thrust::minimum<float>());
    float ymax = thrust::reduce(thrust::device, din_y, din_y + count, -1e8f, thrust::maximum<float>());
    min = std::min(min, ymin);
    max = std::max(max, ymax);
    float zmin = thrust::reduce(thrust::device, din_z, din_z + count, 1e8f, thrust::minimum<float>());
    float zmax = thrust::reduce(thrust::device, din_z, din_z + count, -1e8f, thrust::maximum<float>());
    min = std::min(min, zmin);
    max = std::max(max, zmax);
    int threads = 1024;
    int blocks = (int)ceil((1.0 * count) / threads);
    scale_points_to_unity <<< blocks, threads >> > (din_x, din_y, din_z, min, max, count);
    cudaDeviceSynchronize();

	int times = count / max_threads;
    int start = 0;
    for (int i = 0; i < times; i++) {
        int threads = 1024;
        int blocks = (int)ceil((1.0 * max_threads) / threads);
        calculateEGG<<<blocks, threads>>>(din_x, din_y, din_z, dneighbors, ratio, count, maxNeighbors, start, max_threads);
        cudaDeviceSynchronize();
        start += max_threads;
    }
    if (start < count){
        int threads = 1024;
        int blocks = (int)ceil((1.0 * (count - start)) / threads);
        calculateEGG<<<blocks, threads>>>(din_x, din_y, din_z, dneighbors, ratio, count, maxNeighbors, start, count - start);
        cudaDeviceSynchronize();
    }

    for (int k = 0; k < iterationCount; k++) {
        start = 0;
        for (int i = 0; i < times; i++) {
            int threads = 1024;
            int blocks = (int)ceil((1.0 * (max_threads)) / threads);
            taubin_step<<<blocks, threads>>>(din_x, din_y, din_z,
                                             dout_x, dout_y, dout_z, count, lambda, dneighbors,
                                             maxNeighbors, isRegularized, start, max_threads);
            cudaDeviceSynchronize();
            start += max_threads;
        }
        if (start < count){
            int threads = 1024;
            int blocks = (int)ceil((1.0 * (count - start)) / threads);
            taubin_step<<<blocks, threads>>>(din_x, din_y, din_z,
                                             dout_x, dout_y, dout_z, count, lambda, dneighbors,
                                             maxNeighbors, isRegularized, start, count - start);
            cudaDeviceSynchronize();
        }
        std::swap(din_x, dout_x);
        std::swap(din_y, dout_y);
        std::swap(din_z, dout_z);

        start = 0;
        for (int i = 0; i < times; i++) {
            int threads = 1024;
            int blocks = (int)ceil((1.0 * (max_threads)) / threads);
            taubin_step<<<blocks, threads>>>(din_x, din_y, din_z,
                                             dout_x, dout_y, dout_z, count, mu, dneighbors,
                                             maxNeighbors, isRegularized, start, max_threads);
            cudaDeviceSynchronize();
            start += max_threads;
        }
        if (start < count){
            int threads = 1024;
            int blocks = (int)ceil((1.0 * (count - start)) / threads);
            taubin_step<<<blocks, threads>>>(din_x, din_y, din_z,
                                             dout_x, dout_y, dout_z, count, mu, dneighbors,
                                             maxNeighbors, isRegularized, start, count - start);
            cudaDeviceSynchronize();
        }
        std::swap(din_x, dout_x);
        std::swap(din_y, dout_y);
        std::swap(din_z, dout_z);
    }

    threads = 1024;
	blocks = (int)ceil((1.0 * count) / threads);
    produce_output << <blocks, threads >> > (din_x, din_y, din_z, dout_x, dout_y, dout_z, min, max, count);
    cudaDeviceSynchronize();

    cudaMemcpy(out_x, dout_x, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_y, dout_y, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_z, dout_z, size, cudaMemcpyDeviceToHost);
    cudaFree(din_x); cudaFree(din_y); cudaFree(din_z);
    cudaFree(dout_x); cudaFree(dout_y); cudaFree(dout_z);
    cudaFree(dneighbors);
}
