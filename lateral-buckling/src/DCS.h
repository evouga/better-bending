#ifndef DCS_H
#define DCS_H

#include "Model.h"
#include "../include/MeshConnectivity.h"
#include "../include/MidedgeAngleSinFormulation.h"
#include "../include/ElasticShell.h"
#include "../include/RestState.h"
#include "../include/StVKMaterial.h"

class DCS : public Model
{
public:
	DCS(const Eigen::MatrixXd& restV, const Eigen::MatrixXi& restF,
		double h, double lameAlpha, double lameBeta)
		: mesh_(restF), mat_()
	{
		LibShell::MidedgeAngleSinFormulation::initializeExtraDOFs(restEdgeDOFs_, mesh_, restV);

		// set uniform thicknesses
		restState_.thicknesses.resize(mesh_.nFaces(), h);
		restState_.lameAlpha.resize(mesh_.nFaces(), lameAlpha);
		restState_.lameBeta.resize(mesh_.nFaces(), lameBeta);

		// initialize first and second fundamental forms to those of input mesh
		LibShell::ElasticShell<LibShell::MidedgeAngleSinFormulation>::
			firstFundamentalForms(mesh_, restV, restState_.abars);

		restState_.bbars.resize(mesh_.nFaces());
		for (int i = 0; i < mesh_.nFaces(); i++)
		{
			restState_.bbars[i].setZero();
		}
	}

	virtual const Eigen::VectorXd& restEdgeDOFs() const
	{
		return restEdgeDOFs_;
	}

	virtual double energy(const Eigen::MatrixXd& curPos, const Eigen::VectorXd& curEdgeDOFs, Eigen::VectorXd* deriv, Eigen::SparseMatrix<double>* hess) const
	{
		std::vector<Eigen::Triplet<double> > Hcoeffs;
		double energy = LibShell::ElasticShell<LibShell::MidedgeAngleSinFormulation>::elasticEnergy(mesh_, curPos, curEdgeDOFs, mat_, restState_, deriv, hess ? &Hcoeffs : nullptr, LibShell::HessianProjectType::kNone);
		if (hess)
			hess->setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());
		return energy;
	}

	virtual bool dofsValid(const Eigen::MatrixXd& curPos, const Eigen::VectorXd& curEdgeDOFs) const
	{
		return true;
	}

private:
	LibShell::MeshConnectivity mesh_;
	Eigen::VectorXd restEdgeDOFs_;
	LibShell::StVKMaterial<LibShell::MidedgeAngleSinFormulation> mat_;
	LibShell::MonolayerRestState restState_;
};

#endif