#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <Eigen/Core>
#include <string>

namespace LibShell
{
	class MeshConnectivity;
}

void visualizeEdgeDirectors(const Eigen::MatrixXd& V, const LibShell::MeshConnectivity& mesh, const Eigen::VectorXd& edgeDOFs,
	Eigen::MatrixXd &edgeMidpoints,
	Eigen::MatrixXd &edgeNormals);

void writePointCloud(const std::string& filename, const Eigen::MatrixXd& V, const Eigen::MatrixXd& N);
void readPointCloud(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXd& N);
void writeQuadElements(const std::string& filename, const Eigen::MatrixXd& vertPos, const Eigen::MatrixXd& edgePos, const Eigen::MatrixXd& visV, const Eigen::MatrixXi& visF, const Eigen::VectorXd& visEnergyDensity, const Eigen::MatrixXd& edgeQuadPoints, const Eigen::MatrixXd& edgeNormals);
void readQuadElements(const std::string& filename, Eigen::MatrixXd& elementNodes, Eigen::MatrixXd& visV, Eigen::MatrixXi& visF, Eigen::VectorXd& visEnergyDensity, Eigen::MatrixXd& edgeQuadPoints, Eigen::MatrixXd& edgeNormals);

#endif