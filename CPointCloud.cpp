#include <cstdio>
#include <cstdlib>
#include <map>
#include <iostream>
#include <numeric>

#include "CPointCloud.h"

void CPointCloud::createUmbrellaFromTriangulation(std::vector<std::map<int, UmbrellaElement>>& umbrellas) const
{
	FILE* file = fopen("vase_sor_smoothed.off", "r");
	char off[10];
	int numbPnts, numbTris, numbEdges;
	fscanf(file, "%s", off);
	fscanf(file, "%d %d %d", &numbPnts, &numbTris, &numbEdges);
	//FILE* file1 = fopen("info.txt", "w");
	//fprintf(file1, "%d %d\n", numbPnts, numbTris);
	for (int i = 0; i < numbPnts; i++)
	{
		float x, y, z;
		fscanf(file, "%f %f %f", &x, &y, &z);
	}
	umbrellas.resize(numbPnts);



	for (int t = 0; t < numbTris; t++)
	{
		int numb, i, j, k;
		fscanf(file, "%d %d %d %d", &numb, &i, &j, &k);
		for (int u = 0; u < 3; u++) {
			if (umbrellas[i].find(j) != umbrellas[i].end())
			{
				umbrellas[i][j].indexNext = k;
			}
			else
			{
				umbrellas[i][j] = { j, -1, k };
			}
			if (umbrellas[j].find(i) != umbrellas[j].end())
			{
				umbrellas[j][i].indexPrevious = k;
			}
			else
			{
				umbrellas[j][i] = { i, k, -1 };
			}
			int i1 = i, j1 = j, k1 = k;
			i = j1; j = k1; k = i1;
		}
	}
	fclose(file);
}
spPointCloud CPointCloud::CreatePointCloudFromFileWithNormals(const char* fileName, bool onlyPnts)
{
	int countPoints = 0;
	FILE* file = fopen(fileName, "r");
	CPointEx3D pnt, nrm;
	spPointCloud pc;
	if (!onlyPnts)
	{
		while (fscanf(file, "%lf %lf %lf %lf %lf %lf", &pnt[0], &pnt[1], &pnt[2], &nrm[0], &nrm[1], &nrm[2]) == 6)
			countPoints++;
		rewind(file);
		pc = std::make_shared<CPointCloud>(countPoints);
		int currentIndex = 0;
		while (fscanf(file, "%lf %lf %lf %lf %lf %lf", &pnt[0], &pnt[1], &pnt[2], &nrm[0], &nrm[1], &nrm[2]) == 6)
		{
			(*pc)[currentIndex].first = pnt;
			(*pc)[currentIndex++].second = nrm;
		}
	}
	else
	{
		while (fscanf(file, "%lf %lf %lf", &pnt[0], &pnt[1], &pnt[2]) == 3)
			countPoints++;
		rewind(file);
		pc = std::make_shared<CPointCloud>(countPoints);
		int currentIndex = 0;
		while (fscanf(file, "%lf %lf %lf", &pnt[0], &pnt[1], &pnt[2]) == 3)
		{
			(*pc)[currentIndex].first = pnt;
			(*pc)[currentIndex++].second = nrm;
		}
	}
	fclose(file);
	return pc;
}

CPointCloud::CPointCloud(const char* file, bool onlyPnts)
{
	auto pc = CreatePointCloudFromFileWithNormals(file, onlyPnts);
	*this = *pc;
}

void CPointCloud::savePointCloud(const char* fileName, bool withNormals) const
{
	FILE* file = fopen(fileName, "w");
	for (int i = 0; i < this->GetSize(); i++)
	{
		if (withNormals)
			fprintf(file, "%lf %lf %lf %lf %lf %lf\n", (*this)[i].first[0], (*this)[i].first[1], (*this)[i].first[2],
				normals[i][0], normals[i][1], normals[i][2]);
		else fprintf(file, "%lf %lf %lf\n", (*this)[i].first[0], (*this)[i].first[1], (*this)[i].first[2]);
	}
	fclose(file);
}

CPointCloud::CPointCloud(const CPointCloud& pcSource):
	CVector<CPointEx3D>(pcSource)
{
	normals = pcSource.normals;
}

CPointCloud& CPointCloud::operator=(const CPointCloud& d)
{
	if (this != &d) {
		CVector<CPointEx3D>::operator=(d);
		normals = d.normals;
	}
	return *this;
}

CPointCloud::CPointCloud(int MaxSize, int GrowRate):
	CVector<CPointEx3D>(MaxSize, GrowRate),
	normals(MaxSize, GrowRate)
{}

CPointEx3D CPointCloud::FindBaryCenter() const
{
	CPointEx3D baryCenter(0., 0., 0.);
	for (int i = 0; i < GetSize(); i++)
		baryCenter += (*this)[i].first;
	if (this->GetSize())
		baryCenter /= this->GetSize();
	return baryCenter;
}

CPointEx3D CPointCloud::FindBaryCenter(const CVector<CPointEx3D>& pnts)
{
	CPointEx3D baryCenter(0., 0., 0.);
	for (int i = 0; i < pnts.GetSize(); i++) baryCenter += pnts[i];
	if (pnts.GetSize())
		baryCenter /= pnts.GetSize();
	return baryCenter;
}

void CPointCloud::Add(const CPointEx3D& point)
{
	CVector<CPointEx3D>::Add(point);
	normals.Add(CPointEx3D(0., 0., 0.));
}

void CPointCloud::Add(const CPointEx3D& point, const CPointEx3D& normal)
{
	CVector<CPointEx3D>::Add(point);
	normals.Add(normal);
}

std::pair<CPointEx3D&, CPointEx3D&> CPointCloud::operator[](int i)
{
	CPointEx3D& a = CVector<CPointEx3D>::operator[](i);
	CPointEx3D& b = normals[i];
	return { a, b };
}

std::pair<CPointEx3D, CPointEx3D> CPointCloud::operator[](int i) const
{
	const CPointEx3D& a = CVector<CPointEx3D>::operator[](i);
	const CPointEx3D& b = normals[i];
	return { a, b };
}

double CPointCloud::AreaOfTriangle(const CPointEx3D& pnt1, const CPointEx3D& pnt2, const CPointEx3D& pnt3)
{
	const double a = (pnt2 - pnt1).norm();
	const double b = (pnt3 - pnt2).norm();
	const double c = (pnt3 - pnt1).norm();
	const double s = (a + b + c) / 2;
	return sqrt(s * (s - a) * (s - b) * (s - c));
}

CPointEx3D CPointCloud::NormalOfTriangle(const CPointEx3D& pnt1, const CPointEx3D& pnt2, const CPointEx3D& pnt3)
{
	CPointEx3D vec12 = pnt2 - pnt1;
	CPointEx3D vec13 = pnt3 - pnt1;
	CPointEx3D normal = cross(vec12, vec13);
	if (normal.norm() == 0.0) return CPointEx3D(0, 0, 0);
	normal.normalize();
	return normal;
}

CDoubleMatrix CPointCloud::CreateLocalCoordinateSystem(const CPointEx3D& o, const CPointEx3D& d, Axis axis)
{
	char indexSet[3];
	switch (axis)
	{
		case X_AXIS:
		{
			indexSet[0] = 1;
			indexSet[1] = 2;
			indexSet[2] = 0;
			break;
		}
		case Y_AXIS:
		{
			indexSet[0] = 2;
			indexSet[1] = 0;
			indexSet[2] = 1;
			break;
		}
		case Z_AXIS:
		{
			indexSet[0] = 0;
			indexSet[1] = 1;
			indexSet[2] = 2;
		}
	}
	CPointEx3D firstPerpendicular = GetNormalizedPerpendicularVectorToVector(d);
	CPointEx3D secondPerpendicular = cross(d, firstPerpendicular);
	secondPerpendicular.normalize();
	CDoubleMatrix transformationMatrix(4, 4);
	transformationMatrix.Identity();
	transformationMatrix(indexSet[0], 0) = firstPerpendicular[0];
	transformationMatrix(indexSet[0], 1) = firstPerpendicular[1];
	transformationMatrix(indexSet[0], 2) = firstPerpendicular[2];
	transformationMatrix(indexSet[1], 0) = secondPerpendicular[0];
	transformationMatrix(indexSet[1], 1) = secondPerpendicular[1];
	transformationMatrix(indexSet[1], 2) = secondPerpendicular[2];
	transformationMatrix(indexSet[2], 0) = d[0];
	transformationMatrix(indexSet[2], 1) = d[1];
	transformationMatrix(indexSet[2], 2) = d[2];
	auto t = Mult(transformationMatrix, PointToVector(o));
	transformationMatrix(0, 3) = -t[0];
	transformationMatrix(1, 3) = -t[1];
	transformationMatrix(2, 3) = -t[2];
	return transformationMatrix;
}


bool CPointCloud::sortFunction(sortElement element1, sortElement element2)
{
	return element1.value < element2.value;
}

CVector<double> CPointCloud::PointToVector(const CPointEx3D& pnt, bool homogenous)
{
	CVector<double> vec(4);
	vec[0] = pnt[0];
	vec[1] = pnt[1];
	vec[2] = pnt[2];
	if (homogenous) vec[3] = 1.0;
	else vec[3] = 0.0;
	return vec;
}

CPointEx3D CPointCloud::VectorToPoint(const CVector<double>& vec)
{
	CPointEx3D pnt;
	pnt[0] = vec[0];
	pnt[1] = vec[1];
	pnt[2] = vec[2];
	return pnt;
}

spPointCloud CPointCloud::rotatePointCloud(const CMatrix<double>& mat) const
{
	spPointCloud pc = std::make_shared<CPointCloud>(this->GetSize());
	for (int i = 0; i < this->GetSize(); i++)
	{
		auto pnt = Mult(mat, PointToVector((*this)[i].first));
		auto nrm = Mult(mat, PointToVector((*this)[i].second, false));
		(*pc)[i].first = VectorToPoint(pnt);
		(*pc)[i].second = VectorToPoint(nrm);
	}
	return pc;
}

CPointEx3D* CPointCloud::begin()
{
	return m_pData;
}

CPointEx3D* CPointCloud::end()
{
	if (m_pData)
		return &m_pData[this->GetSize()];
	return nullptr;
}

const CPointEx3D* CPointCloud::begin() const
{
	return m_pData;
}

const CPointEx3D* CPointCloud::end() const
{
	if (m_pData)
		return &m_pData[this->GetSize()];
	return nullptr;
}

void CPointCloud::scaleToUnity(CVector<CPointEx3D>& points)
{
	double maxx = -1e10;
	double minx =  1e10;
	double maxy = -1e10;
	double miny = 1e10;
	double maxz = -1e10;
	double minz = 1e10;
	for (int i = 0; i < points.GetSize(); i++)
	{
		if (points[i].x > maxx) maxx = points[i].x;
		if (points[i].x < minx) minx = points[i].x;
		if (points[i].y > maxy) maxy = points[i].y;
		if (points[i].y < miny) miny = points[i].y;
		if (points[i].z > maxz) maxz = points[i].z;
		if (points[i].z < minz) minz = points[i].z;
	}

	double max = maxx;
	if (max < maxy) max = maxy;
	if (max < maxz) max = maxz;

	double min = minx;
	if (min > miny) min = miny;
	if (min > minz) min = minz;
	for (int i = 0; i < points.GetSize(); i++)
	{
		points[i].x = (points[i].x - min) / (max - min);
		points[i].y = (points[i].y - min) / (max - min);
		points[i].z = (points[i].z - min) / (max - min);
	}

	auto bc = FindBaryCenter(points);
	for (int i = 0; i < points.GetSize(); i++)
	{
		points[i] -= bc;
	}
}

CPointEx3D CPointCloud::GetNormalizedPerpendicularVectorToVector(const CPointEx3D& inVector)
{
	double max = fabs(inVector[0]);
	int cordIndex = 0;

	if (max < fabs(inVector[1]))
	{
		cordIndex = 1;
		max = fabs(inVector[1]);
	}

	if (max < fabs(inVector[2]))
	{
		cordIndex = 2;
	}
	CPointEx3D outVector;
	outVector[0] = 1.0;
	outVector[1] = 1.0;
	outVector[2] = 1.0;

	switch (cordIndex)
	{
	case 0:
		outVector[0] = (-inVector[1] * outVector[1] - inVector[2] * outVector[2]) / inVector[0];
		break;
	case 1:
		outVector[1] = (-inVector[0] * outVector[0] - inVector[2] * outVector[2]) / inVector[1];
		break;
	case 2:
		outVector[2] = (-inVector[0] * outVector[0] - inVector[1] * outVector[1]) / inVector[2];
		break;
	}
	outVector.normalize();
	return outVector;
}

CDoubleMatrix CPointCloud::CreatePatchLocalCoordinateSystem(const CPointEx3D& p, const CPointEx3D& n, 
	CPointEx3D& firstPerpendicular, CPointEx3D& secondPerpendicular)
{
	firstPerpendicular = GetNormalizedPerpendicularVectorToVector(n);
	secondPerpendicular = cross(n, firstPerpendicular);
	secondPerpendicular.normalize();
	CDoubleMatrix transformationMatrix(4, 4);
	transformationMatrix.Identity();
	transformationMatrix(0, 0) = firstPerpendicular[0];
	transformationMatrix(0, 1) = firstPerpendicular[1];
	transformationMatrix(0, 2) = firstPerpendicular[2];
	transformationMatrix(1, 0) = secondPerpendicular[0];
	transformationMatrix(1, 1) = secondPerpendicular[1];
	transformationMatrix(1, 2) = secondPerpendicular[2];
	transformationMatrix(2, 0) = n[0];
	transformationMatrix(2, 1) = n[1];
	transformationMatrix(2, 2) = n[2];
	auto t = Mult(transformationMatrix, PointToVector(p));
	transformationMatrix(0, 3) = -t[0];
	transformationMatrix(1, 3) = -t[1];
	transformationMatrix(2, 3) = -t[2];
	return transformationMatrix;
}
