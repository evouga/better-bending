#ifndef S1SIN_H
#define S1SIN_H

#include "Model.h"
#include "../include/MidedgeAngleSinFormulation.h"
#include "../include/RestState.h"
#include "../include/ElasticShell.h"
#include "../include/StVKMaterial.h"

class S1Sin : public Model
{
private:
    void optimizeEdgeDOFs(const LibShell::MeshConnectivity& mesh, 
        const Eigen::MatrixXd& curPos, 
        Eigen::VectorXd& edgeDOFs,
        const LibShell::StVKMaterial<LibShell::MidedgeAngleSinFormulation> &mat,
        const LibShell::MonolayerRestState &restState
        )
    {
        double tol = 1e-5;
        int nposdofs = curPos.rows() * 3;
        int nedgedofs = edgeDOFs.size();

        std::vector<Eigen::Triplet<double> > Pcoeffs;
        std::vector<Eigen::Triplet<double> > Icoeffs;
        for (int i = 0; i < nedgedofs; i++)
        {
            Pcoeffs.push_back({ i, nposdofs + i, 1.0 });
            Icoeffs.push_back({ i, i, 1.0 });
        }
        Eigen::SparseMatrix<double> P(nedgedofs, nposdofs + nedgedofs);
        P.setFromTriplets(Pcoeffs.begin(), Pcoeffs.end());
        Eigen::SparseMatrix<double> I(nedgedofs, nedgedofs);
        I.setFromTriplets(Icoeffs.begin(), Icoeffs.end());

        Eigen::SparseMatrix<double> PT = P.transpose();

        double reg = 1e-6;
        while (true)
        {
            std::vector<Eigen::Triplet<double> > Hcoeffs;
            Eigen::VectorXd F;
            double origEnergy = LibShell::ElasticShell<LibShell::MidedgeAngleSinFormulation>::elasticEnergy(mesh, curPos, edgeDOFs, mat, restState, &F, &Hcoeffs, LibShell::HessianProjectType::kNone);
            Eigen::VectorXd PF = P * F;
            std::cout << "Force resid now: " << PF.norm() << ", energy: " << origEnergy << ", reg: " << reg;
            if (PF.norm() < tol)
            {
                std::cout << "; converged" << std::endl;
                return;
            }

            std::vector<Eigen::Triplet<double> > Hfiltered;
            for (auto& it : Hcoeffs)
            {
                if (it.row() >= nposdofs && it.col() >= nposdofs)
                {
                    Hfiltered.push_back({ it.row() - nposdofs, it.col() - nposdofs, it.value() });
                }
            }
            Hcoeffs.clear();
            Eigen::SparseMatrix<double> PHPT(nedgedofs, nedgedofs);
            PHPT.setFromTriplets(Hfiltered.begin(), Hfiltered.end());

            Eigen::SparseMatrix<double> M = PHPT + reg * I;
            Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver(M);
            Eigen::VectorXd update = solver.solve(-PF);
            if (solver.info() != Eigen::Success) {
                std::cout << "; solve failed" << std::endl;
                reg *= 2.0;
                continue;
            }
            else
            {
                std::cout << ", Newton decrement: " << update.norm() << std::endl;
            }
            Eigen::VectorXd newedgeDOFs = edgeDOFs + update;
            double newenergy = LibShell::ElasticShell<LibShell::MidedgeAngleSinFormulation>::elasticEnergy(mesh, curPos, newedgeDOFs, mat, restState, NULL, NULL);
            if (newenergy > origEnergy)
            {
                std::cout << "Not a descent step, " << origEnergy << " -> " << newenergy << std::endl;
                reg *= 2.0;
                continue;
            }
            edgeDOFs = newedgeDOFs;
            reg *= 0.5;
        }
    }

    double elasticEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd *restPos,
        const Eigen::MatrixXd& curPos,
        double thickness,
        double lameAlpha,
        double lameBeta)
    {
        Eigen::VectorXd edgeDOFs;
		LibShell::MidedgeAngleSinFormulation::initializeExtraDOFs(edgeDOFs, mesh, restPos ? *restPos : curPos);

        // initialize the rest geometry of the shell
        LibShell::MonolayerRestState restState;

        // set uniform thicknesses
        restState.thicknesses.resize(mesh.nFaces(), thickness);
        restState.lameAlpha.resize(mesh.nFaces(), lameAlpha);
        restState.lameBeta.resize(mesh.nFaces(), lameBeta);

        LibShell::StVKMaterial<LibShell::MidedgeAngleSinFormulation> mat;

        // initialize first and second fundamental forms to those of input mesh
        LibShell::ElasticShell<LibShell::MidedgeAngleSinFormulation>::
            firstFundamentalForms(mesh, restPos ? *restPos : curPos, restState.abars);
        restState.bbars.resize(mesh.nFaces());
        for (int i = 0; i < mesh.nFaces(); i++)
        {
            restState.bbars[i].setZero();
        }

		// if there is a rest state, solve for the edge directors on the rest state
        if (restPos)
        {
            Eigen::VectorXd restEdgeDOFs = edgeDOFs;
            optimizeEdgeDOFs(mesh, *restPos, restEdgeDOFs, mat, restState);
            LibShell::ElasticShell<LibShell::MidedgeAngleSinFormulation>::
                secondFundamentalForms(mesh, *restPos, restEdgeDOFs, restState.bbars);
        }
        
        optimizeEdgeDOFs(mesh, curPos, edgeDOFs, mat, restState);

        return LibShell::ElasticShell<LibShell::MidedgeAngleSinFormulation>::elasticEnergy(mesh, curPos, edgeDOFs, mat, restState, NULL, NULL);
    }

public:
    virtual double measureHalfCylinderEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& restPos,
        const Eigen::MatrixXd& restEdgePos,
        const Eigen::MatrixXd& curPos,
        const Eigen::MatrixXd& curEdgePos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        double radius,
        double height)
    {
		return elasticEnergy(mesh, &restPos, curPos, thickness, lameAlpha, lameBeta);
    }

    virtual double measureRestCurvedHalfCylinderEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& restPos,
        const Eigen::MatrixXd& restEdgePos,
        const Eigen::MatrixXd& curPos,
        const Eigen::MatrixXd& curEdgePos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        double radius,
        double height)
    {
        return elasticEnergy(mesh, &restPos, curPos, thickness, lameAlpha, lameBeta);
    }

    virtual double measureVeeEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& restPos,
        const Eigen::MatrixXd& restEdgePos,
        const Eigen::MatrixXd& curPos,
        const Eigen::MatrixXd& curEdgePos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        double radius,
        double height)
    {
        return elasticEnergy(mesh, &restPos, curPos, thickness, lameAlpha, lameBeta);
    }

    virtual double measureSphereEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& curPos,
        const Eigen::MatrixXd& curEdgePos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        double radius)
    {
        return elasticEnergy(mesh, NULL, curPos, thickness, lameAlpha, lameBeta);
    }

    virtual double measureHalfTorusEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& curPos,
        const Eigen::MatrixXd& curEdgePos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        double innerRadius,
        double outerRadius)
    {
        return elasticEnergy(mesh, NULL, curPos, thickness, lameAlpha, lameBeta);
	}

    virtual double measureHelixEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& restPos, 
        const Eigen::MatrixXd& restEdgePos,
        const Eigen::MatrixXd& curPos,
        const Eigen::MatrixXd& curEdgePos,
        double thickness,
        double lameAlpha,
        double lameBeta,
		double helixRadius,
		double helixHeight,
        double helixSurfaceHeight)
    {
        return elasticEnergy(mesh, &restPos, curPos, thickness, lameAlpha, lameBeta);
	}
};

#endif
