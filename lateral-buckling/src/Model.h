#ifndef MODEL_H
#define MODEL_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class Model
{
public:
	virtual const Eigen::VectorXd& restEdgeDOFs() const = 0;

	virtual double energy(const Eigen::MatrixXd& curPos, const Eigen::VectorXd& curEdgeDOFs, Eigen::VectorXd *deriv, Eigen::SparseMatrix<double> *hess) const = 0;

	virtual bool dofsValid(const Eigen::MatrixXd& curPos, const Eigen::VectorXd& curEdgeDOFs) const = 0;
};

#endif