#include "Visualization.h"
#include "../include/MeshConnectivity.h"
#include "../src/GeometryDerivatives.h"
#include <fstream>
#include <cassert>
#include <Eigen/Dense>

void visualizeEdgeDirectors(const Eigen::MatrixXd& V, const LibShell::MeshConnectivity& mesh, const Eigen::VectorXd& edgeDOFs,
	Eigen::MatrixXd &edgeMidpoints,
	Eigen::MatrixXd &edgeNormals)
{
	int nEdges = mesh.nEdges();
	assert(edgeDOFs.size() == nEdges);
	edgeMidpoints.resize(nEdges, 3);
	edgeNormals.resize(nEdges, 3);
	for (int e = 0; e < nEdges; e++)
	{
		int v0 = mesh.edgeVertex(e, 0);
		int v1 = mesh.edgeVertex(e, 1);
		Eigen::Vector3d p0 = V.row(v0).transpose();
		Eigen::Vector3d p1 = V.row(v1).transpose();
		Eigen::Vector3d evec = (p1 - p0).normalized();
		Eigen::Vector3d midpoint = 0.5 * (p0 + p1);

		double theta = 0;
		int f;
		
		if (mesh.edgeFace(e, 1) != -1 && mesh.edgeFace(e, 0) != -1)
		{
			int f0 = mesh.edgeFace(e, 0);
			Eigen::Vector3d fNormal0 = LibShell::faceNormal(mesh, V, f0, 0, NULL, NULL).normalized();
			int f1 = mesh.edgeFace(e, 1);
			Eigen::Vector3d fNormal1 = LibShell::faceNormal(mesh, V, f1, 0, NULL, NULL).normalized();
			theta = LibShell::angle(fNormal0, fNormal1, evec, NULL, NULL);
			f = f0;			
		}
		else if(mesh.edgeFace(e, 0) != -1)
		{
			f = mesh.edgeFace(e, 0);
		}
		else if (mesh.edgeFace(e, 1) != -1)
		{
			f = mesh.edgeFace(e, 1);
		}
		else
		{
			exit(-1); // isolated edge?
		}

		Eigen::Vector3d fNormal = LibShell::faceNormal(mesh, V, f, 0, NULL, NULL).normalized();
		
		int edgeidx = 0;
		for (int i = 0; i < 3; i++)
		{
			if (mesh.faceEdge(f, i) == e)
			{
				edgeidx = i;
			}
		}
		double orient = mesh.faceEdgeOrientation(f, edgeidx) == 0 ? 1.0 : -1.0;
		double alpha = 0.5 * theta + orient * edgeDOFs[e];

		Eigen::Vector3d director = -std::sin(alpha) * orient * fNormal.cross(evec) + std::cos(alpha) * fNormal;
		edgeMidpoints.row(e) = midpoint;
		edgeNormals.row(e) = director.transpose();
	}
}

void writePointCloud(const std::string& filename, const Eigen::MatrixXd& V, const Eigen::MatrixXd& N)
{
	std::ofstream ofs(filename);
	int nverts = V.rows();
	assert(V.cols() == 3);
	assert(N.rows() == nverts);
	assert(N.cols() == 3);
	ofs << nverts << std::endl;
	for(int i=0; i < nverts; i++)
	{
		ofs << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << " "
			<< N(i, 0) << " " << N(i, 1) << " " << N(i, 2) << std::endl;
	}
}

void readPointCloud(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXd& N)
{
	std::ifstream ifs(filename);
	int nverts;
	ifs >> nverts;
	assert(ifs);
	assert(nverts >= 0);
	V.resize(nverts, 3);
	N.resize(nverts, 3);
	for(int i=0; i< nverts; i++)
	{
		ifs >> V(i, 0) >> V(i, 1) >> V(i, 2)
			>> N(i, 0) >> N(i, 1) >> N(i, 2);
	}
	assert(ifs);
	std::string dummy;
	ifs >> dummy;
	assert(!ifs);
}

void writeQuadElements(const std::string& filename, const Eigen::MatrixXd& vertPos, const Eigen::MatrixXd& edgePos, 
	const Eigen::MatrixXd& visV, const Eigen::MatrixXi& visF, const Eigen::VectorXd& visEnergyDensity,
	const Eigen::MatrixXd& edgeQuadPoints, const Eigen::MatrixXd& edgeNormals)
{
	std::ofstream ofs(filename);
	int nverts = vertPos.rows() + edgePos.rows();
	ofs << nverts << std::endl;
	for(int i=0; i<vertPos.rows(); i++)
		ofs << vertPos(i, 0) << " " << vertPos(i, 1) << " " << vertPos(i, 2) << std::endl;
	for(int i=0; i<edgePos.rows(); i++)
		ofs << edgePos(i, 0) << " " << edgePos(i, 1) << " " << edgePos(i, 2) << std::endl;
	int nvverts = visV.rows();
	int nvfaces = visF.rows();
	ofs << nvverts << " " << nvfaces << std::endl;
	for(int i=0; i<nvverts; i++)
		ofs << visV(i, 0) << " " << visV(i, 1) << " " << visV(i, 2) << " " << visEnergyDensity[i] << std::endl;
	for (int i = 0; i < nvfaces; i++)
		ofs << visF(i, 0) << " " << visF(i, 1) << " " << visF(i, 2) << std::endl;
	int nqpts = edgeQuadPoints.rows();
	ofs << nqpts << std::endl;
	for (int i = 0; i < nqpts; i++)
	{
		ofs << edgeQuadPoints(i, 0) << " " << edgeQuadPoints(i, 1) << " " << edgeQuadPoints(i, 2) << " " << edgeNormals(i, 0) << " " << edgeNormals(i, 1) << " " << edgeNormals(i, 2) << std::endl;
	}
}

void readQuadElements(const std::string &filename, Eigen::MatrixXd &elementNodes,
	Eigen::MatrixXd &visV, Eigen::MatrixXi &visF, Eigen::VectorXd &visEnergyDensity,
	Eigen::MatrixXd &edgeQuadPoints, Eigen::MatrixXd &edgeNormals)
{
	std::ifstream ifs(filename);
	int nverts;
	ifs >> nverts;
	elementNodes.resize(nverts, 3);
	for (int i = 0; i < nverts; i++)
	{
		ifs >> elementNodes(i, 0) >> elementNodes(i, 1) >> elementNodes(i, 2);
	}
	int nvverts, nvfaces;
	ifs >> nvverts >> nvfaces;
	visV.resize(nvverts, 3);
	visF.resize(nvfaces, 3);
	visEnergyDensity.resize(nvverts);
	for (int i = 0; i < nvverts; i++)
	{
		ifs >> visV(i, 0) >> visV(i, 1) >> visV(i, 2) >> visEnergyDensity[i];
	}
	for (int i = 0; i < nvfaces; i++)
	{
		ifs >> visF(i, 0) >> visF(i, 1) >> visF(i, 2);
	}	
	int neqpts;
	ifs >> neqpts;
	edgeQuadPoints.resize(neqpts, 3);
	edgeNormals.resize(neqpts, 3);
	for (int i = 0; i < neqpts; i++)
	{
		ifs >> edgeQuadPoints(i, 0) >> edgeQuadPoints(i, 1) >> edgeQuadPoints(i, 2) >> edgeNormals(i, 0) >> edgeNormals(i, 1) >> edgeNormals(i, 2);
	}
}