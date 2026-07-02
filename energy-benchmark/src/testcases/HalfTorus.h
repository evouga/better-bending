#ifndef HALFTORUS_H
#define HALFTORUS_H

#include <Eigen/Core>
#include <vector>

namespace LibShell
{
    class MeshConnectivity;
};

void makeHalfTorus(bool regular, double innerRadius, double outerRadius, double triangleArea,
    Eigen::MatrixXd& V,
    Eigen::MatrixXd &edgeV,
    LibShell::MeshConnectivity &mesh);


#endif