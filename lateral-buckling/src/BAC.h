#ifndef BAC_H
#define BAC_H

#include "Model.h"
#include "../include/MeshConnectivity.h"
#include "../include/MidedgeAngleTanFormulation.h"
#include "../include/ElasticShell.h"
#include "../include/RestState.h"
#include "../include/StVKMaterial.h"

class BAC : public Model
{
public:
	BAC(const Eigen::MatrixXd &restV, const Eigen::MatrixXi& restF,
		double h, double lameAlpha, double lameBeta) 
		: mesh_(restF), mat_()
	{
		LibShell::MidedgeAngleTanFormulation::initializeExtraDOFs(restEdgeDOFs_, mesh_, restV);

		// set uniform thicknesses
		restState_.thicknesses.resize(mesh_.nFaces(), h);
		restState_.lameAlpha.resize(mesh_.nFaces(), lameAlpha);
		restState_.lameBeta.resize(mesh_.nFaces(), lameBeta);

		// initialize first and second fundamental forms to those of input mesh
		LibShell::ElasticShell<LibShell::MidedgeAngleTanFormulation>::
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
		double energy = LibShell::ElasticShell<LibShell::MidedgeAngleTanFormulation>::elasticEnergy(mesh_, curPos, curEdgeDOFs, mat_, restState_, deriv, hess ? &Hcoeffs : nullptr, LibShell::HessianProjectType::kNone);
		if(hess)
			hess->setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());
		return energy;
	}

	virtual bool dofsValid(const Eigen::MatrixXd& curPos, const Eigen::VectorXd& curEdgeDOFs) const
	{
		return LibShell::MidedgeAngleTanFormulation::edgeDOFsValid(mesh_, curPos, curEdgeDOFs);
	}

private:	
	LibShell::MeshConnectivity mesh_;
	Eigen::VectorXd restEdgeDOFs_;
	LibShell::StVKMaterial<LibShell::MidedgeAngleTanFormulation> mat_;
	LibShell::MonolayerRestState restState_;
};

#endif