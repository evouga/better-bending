#ifndef VEE_H
#define VEE_H

#include <Eigen/Core>
#include <vector>

namespace LibShell
{
    class MeshConnectivity;
};

void makeVee(double width, double height, double creaseAngle, double triangleArea,
    Eigen::MatrixXd& flatV,
    Eigen::MatrixXd& flatEdgeV,
    Eigen::MatrixXd& V,
    Eigen::MatrixXd& edgeV,
    LibShell::MeshConnectivity& mesh);

#endif