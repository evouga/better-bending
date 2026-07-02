#ifndef SPHERE_H
#define SPHERE_H

#include <Eigen/Core>

namespace LibShell
{
	class MeshConnectivity;
};

void makeSphere(double radius, double triangleArea, Eigen::MatrixXd& V, Eigen::MatrixXd &edgeV, LibShell::MeshConnectivity &mesh);
void sphereFromSamples(int samples, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

#endif
