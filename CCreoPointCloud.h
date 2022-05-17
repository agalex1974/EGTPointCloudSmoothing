#ifndef CREO_POINT_CLOUD_H
#define CREO_POINT_CLOUD_H
#include <queue>
#include "CPointCloud.h"
// CGAL HEADERS///////////////////////////////////
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <boost/iterator/zip_iterator.hpp>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/property_map.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/edge_aware_upsample_point_set.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
///////////////////////////////////////////////
/// <summary>
// CGAL DEFINITIONS
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef boost::tuple<Point_3, int> Point_and_int;
typedef CGAL::Search_traits_3<Kernel>                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
	CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
	Traits_base> Traits;
typedef CGAL::Random_points_in_cube_3<Point_3>       Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;
typedef CGAL::Kd_tree<Traits> Tree;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;
typedef CGAL::Fuzzy_iso_box<Traits> Fuzzy_iso_box;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
typedef K_neighbor_search::Distance                         Distance;
typedef CGAL::Parallel_if_available_tag Concurrency_tag;
// Point with normal vector stored in a std::pair.
typedef std::pair<Point_3, Vector_3> CGALPointVectorPair;
typedef std::vector<CGALPointVectorPair> CGALPointList;

/// </summary>

using PrincipalCurvatures = std::pair<std::shared_ptr<std::vector<double>>, std::shared_ptr<std::vector<double>>>;

class CCreoPointCloud : public CPointCloud
{
public:
	explicit CCreoPointCloud(const char* pntCloudFile, bool onlyPnts = false);
	~CCreoPointCloud();
	std::vector<int> GetKNearestNeighborIndex(const CPointEx3D& pnt, int K) const;
	bool isEllipicGabrielNeighbor(int i, const std::vector<int>& NNs, double a) const;
	std::vector<int> GabrielkNeighborhood(int i, int k, const std::vector<std::vector<int>>& EGG);
    int calculateEGG(std::vector<std::vector<int>>& EGG, double ratio) const;
    void calculateNeighbors(std::vector<int>& NNs, int neighborCount) const;
    void smooth_iterative(int nb_iter, double lambda, double mu, double gabrielRatio, int isRegularized);
private:
	void updateSpatialTree();
	std::shared_ptr<Tree> pointCloudSpatialTree;
};

#endif