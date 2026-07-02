#ifndef HELIX_H
#define HELIX_H

#include <Eigen/Core>
#include <vector>

namespace LibShell
{
    class MeshConnectivity;
};

void makeHelix(bool regular, double helixRadius, double helixHeight, double surfaceHeight, double triangleArea,
    Eigen::MatrixXd& flatV,
    Eigen::MatrixXd &flatEdgeV,
    Eigen::MatrixXd& V,
    Eigen::MatrixXd &edgeV,
    LibShell::MeshConnectivity &mesh);


#endif