#ifndef HALFCYLINDER_H
#define HALFCYLINDER_H

#include <Eigen/Core>
#include <vector>

namespace LibShell
{
    class MeshConnectivity;
};

void makeHalfCylinder(bool regular, bool restCurved, double radius, double height, double aniso, double triangleArea,
    Eigen::MatrixXd& restV,
	Eigen::MatrixXd& restEdgeV,
    Eigen::MatrixXd& V,
	Eigen::MatrixXd& edgeV,
    LibShell::MeshConnectivity& mesh);


void getBoundaries(const Eigen::MatrixXi& F, std::vector<int>& bdryVertices);

#endif