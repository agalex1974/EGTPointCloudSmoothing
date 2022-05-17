#include "CCreoPointCloud.h"
#include <vector>
#include <set>
#include <queue>
//#include <mkl.h>
#include <omp.h>
#include <stack>
#include <chrono>
#include <algorithm>
#include "EGT.cuh"

//a *.xyz (coordinate,normal) file
CCreoPointCloud::CCreoPointCloud(const char* pntCloudFile, bool onlyPnts):
	CPointCloud(pntCloudFile, onlyPnts)
{
	std::vector<Point_3> points(GetSize());
	std::vector<int> indices(GetSize());
#pragma omp parallel for
	for (int i = 0; i < GetSize(); i++)
	{
		points[i] = Point_3((*this)[i].first[0], (*this)[i].first[1], (*this)[i].first[2]);
		indices[i] = i;
	}
	pointCloudSpatialTree = std::make_shared<Tree>(boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())),
		boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));
}

std::vector<int> CCreoPointCloud::GabrielkNeighborhood(int i, int k, const std::vector<std::vector<int>>& EGG)
{
	std::vector<int> output;
	std::set<int> ringsSet;
	std::queue <int> neighborQueue;
	std::map<int, int> neighborMap;

	output.push_back(i);
	ringsSet.insert(i);
	neighborQueue.push(i);
	neighborMap[i] = 0;

	while(!neighborQueue.empty())
	{
		int index = neighborQueue.front();
		neighborQueue.pop();
		if (neighborMap[index] < k) {
			for (const auto& el : EGG[index])
			{
				if (ringsSet.insert(el).second)
				{
					output.push_back(el);
					neighborQueue.push(el);
					neighborMap[el] = neighborMap[index] + 1;
				}
			}
		}
	}	
	return output;
}

//i is the index order in NNs not the index of the queryPoint
//we need to check all points 1..i-1
bool CCreoPointCloud::isEllipicGabrielNeighbor(int i, const std::vector<int>& NNs, double a) const
{
	const CPointEx3D& p  = (*this)[NNs[0]].first;
	const CPointEx3D& qi = (*this)[NNs[i]].first;
	const CPointEx3D origin = 0.5 * (p + qi);
	const double d = (qi - p).norm() / 2.0;
	CPointEx3D localXaxis = qi - p;
	localXaxis.normalize();
	CDoubleMatrix transformationMatrix = CreateLocalCoordinateSystem(origin, localXaxis, X_AXIS);
	for (int j = 1; j < i; j++)
	{
		CPointEx3D pnt = VectorToPoint(Mult(transformationMatrix, PointToVector((*this)[NNs[j]].first)));
		double x = pnt.x, y = pnt.y, z = pnt.z;
		double ellipsoidValue = x * x + y * y / (a * a) + z * z / (a * a);
		if (ellipsoidValue < d * d) return false;
	}
	return true;
}

void CCreoPointCloud::smooth_iterative(int nb_iter, double lambda, double mu, double gabrielRatio, int isRegularized)
{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    updateSpatialTree();
	unsigned nb_vertices = unsigned(GetSize());

    std::vector<float> in_x(nb_vertices);
    std::vector<float> in_y(nb_vertices);
    std::vector<float> in_z(nb_vertices);
    std::vector<float> out_x(nb_vertices);
    std::vector<float> out_y(nb_vertices);
    std::vector<float> out_z(nb_vertices);

#pragma omp parallel for
	for (int i = 0; i < GetSize(); i++)
	{
		in_x[i] = (*this)[i].first.x;
        in_y[i] = (*this)[i].first.y;
        in_z[i] = (*this)[i].first.z;
	}
    int maximum_neighbors = 40;
    std::vector<int> neighbors(maximum_neighbors * nb_vertices);
    calculateNeighbors(neighbors, maximum_neighbors);

    EGTsmoothing(in_x.data(), in_y.data(), in_z.data(), nb_vertices, lambda, mu, neighbors.data(), maximum_neighbors,
                 out_x.data(), out_y.data(), out_z.data(), nb_iter, isRegularized, gabrielRatio);
#pragma omp parallel for
	for (int i = 0; i < GetSize(); i++)
	{
		(*this)[i].first.x = out_x[i];
        (*this)[i].first.y = out_y[i];
        (*this)[i].first.z = out_z[i];
	}
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    FILE* file = fopen("TaubinSmoothingTime.txt", "w");
    fprintf(file, "Time taken:%lf", time);
    fclose(file);
	updateSpatialTree();
}

//We assume that the EGG vector is initialized
int CCreoPointCloud::calculateEGG(std::vector<std::vector<int>>& EGG, double ratio) const
{
	int pointsCount = EGG.size();
    int maximum = 0;
#pragma omp parallel for
	for (int i = 0; i < pointsCount; i++)
	{
		auto NNs = GetKNearestNeighborIndex((*this)[i].first, 40);
		for (int j = 1; j < NNs.size(); j++)
		{
			if (isEllipicGabrielNeighbor(j, NNs, ratio))
				EGG[i].push_back(NNs[j]);
		}
        maximum = max((int)EGG[i].size(), maximum);
	}
    return maximum;
}

void CCreoPointCloud::calculateNeighbors(std::vector<int>& NNs, int neighborCount) const
{
    int pointsCount = GetSize();
#pragma omp parallel for
    for (int i = 0; i < pointsCount; i++)
    {
        auto pointNNs = GetKNearestNeighborIndex((*this)[i].first, neighborCount + 1);
        for (int j = 1; j < pointNNs.size(); j++) {
            NNs[(j - 1) * pointsCount + i] = pointNNs[j];
        }
    }
}

std::vector<int> CCreoPointCloud::GetKNearestNeighborIndex(const CPointEx3D& pnt, int K) const
{
	std::vector<int> indexes;
	const Point_3 point(pnt.x, pnt.y, pnt.z);
	K_neighbor_search search(*pointCloudSpatialTree, point, K);
	for (K_neighbor_search::iterator it = search.begin(); it != search.end(); it++) 
	{
		indexes.push_back(boost::get<1>(it->first));
	}
	return indexes;
}

CCreoPointCloud::~CCreoPointCloud()
{
}

void CCreoPointCloud::updateSpatialTree()
{
	std::vector<Point_3> points(GetSize());
	std::vector<int> indices(GetSize());
#pragma omp parallel for
	for (int i = 0; i < GetSize(); i++)
	{
		points[i] = Point_3((*this)[i].first[0], (*this)[i].first[1], (*this)[i].first[2]);
		indices[i] = i;
	}
	pointCloudSpatialTree = std::make_shared<Tree>(boost::make_zip_iterator(boost::make_tuple(points.begin(), indices.begin())),
		boost::make_zip_iterator(boost::make_tuple(points.end(), indices.end())));
}
