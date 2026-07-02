#ifndef QBCR_H
#define QBCR_H

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "Model.h"
#include "../include/MeshConnectivity.h"
#include "../include/MidedgeAngleTanFormulation.h"
#include "../include/ElasticShell.h"
#include "../include/RestState.h"
#include "../include/StVKMaterial.h"

class QBCR : public Model
{
public:
    QBCR(const Eigen::MatrixXd& restV, const Eigen::MatrixXi& restF,
        double h, double lameAlpha, double lameBeta)
        : mesh_(restF), restPos_(restV), mat_()
    {
		buildQBCRBendingMatrix(restV, h, lameAlpha, lameBeta);

        restEdgeDOFs_ = Eigen::VectorXd(0);

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
        int nposdofs = curPos.rows() * 3;

		std::vector<Eigen::Triplet<double> > Hcoeffs;
        Eigen::VectorXd D;
        Eigen::VectorXd edgeDOFs;
        LibShell::MidedgeAngleTanFormulation::initializeExtraDOFs(edgeDOFs, mesh_, restPos_);
        double energy = LibShell::ElasticShell<LibShell::MidedgeAngleTanFormulation>::elasticEnergy(mesh_, curPos, edgeDOFs, mat_, restState_, LibShell::ElasticShell<LibShell::MidedgeAngleTanFormulation>::ET_STRETCHING, deriv ? &D : nullptr, hess ? &Hcoeffs : nullptr, LibShell::HessianProjectType::kNone);
        
        Eigen::VectorXd q = curPos.reshaped<Eigen::RowMajor>();
        energy += 0.5 * q.transpose() * bendingM_ * q;
        if (deriv)
        {
            *deriv = D.segment(0, nposdofs) + bendingM_ * q;
        }

        if (hess)
        {
            for (int k = 0; k < bendingM_.outerSize(); ++k)
            {
                for (Eigen::SparseMatrix<double>::InnerIterator it(bendingM_, k); it; ++it)
                {
                    Hcoeffs.push_back({ (int)it.row(), (int)it.col(), it.value() });
                }
            }
			hess->setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());
        }

        return energy;
    }

    virtual bool dofsValid(const Eigen::MatrixXd& curPos, const Eigen::VectorXd& curEdgeDOFs) const
    {
        return true;
	}

private:
    LibShell::MeshConnectivity mesh_;
    Eigen::MatrixXd restPos_;
    Eigen::VectorXd restEdgeDOFs_;
    Eigen::SparseMatrix<double> bendingM_;
    LibShell::StVKMaterial<LibShell::MidedgeAngleTanFormulation> mat_;
    LibShell::MonolayerRestState restState_;

    static double cotan(Eigen::Vector3d e0, Eigen::Vector3d e1)
    {
        double dot = e0.dot(e1);
        double crossNorm = e0.cross(e1).norm();
        if (crossNorm == 0.0)
            return 0.0;
        return dot / crossNorm;
    }

    void buildQBCRBendingMatrix(       
        const Eigen::MatrixXd& restPos,
        double thickness,
        double lameAlpha,
        double lameBeta)
    {
        int nverts = restPos.rows();
        int nfaces = mesh_.nFaces();
        int nedges = mesh_.nEdges();

        double weight = thickness * thickness * thickness / 12.0 * (lameAlpha + 2.0 * lameBeta);
        std::vector<Eigen::Triplet<double>> Qentries;
        for (int i = 0; i < nedges; i++)
        {
            if (mesh_.edgeOppositeVertex(i, 0) == -1 || mesh_.edgeOppositeVertex(i, 1) == -1)
                continue; // skip boundary edges    
            int v[4];
            v[0] = mesh_.edgeVertex(i, 0);
            v[1] = mesh_.edgeVertex(i, 1);
            v[2] = mesh_.edgeOppositeVertex(i, 0);
            v[3] = mesh_.edgeOppositeVertex(i, 1);
            Eigen::Vector3d e0 = restPos.row(v[1]) - restPos.row(v[0]);
            Eigen::Vector3d e1 = restPos.row(v[2]) - restPos.row(v[0]);
            Eigen::Vector3d e2 = restPos.row(v[3]) - restPos.row(v[0]);
            Eigen::Vector3d e3 = restPos.row(v[2]) - restPos.row(v[1]);
            Eigen::Vector3d e4 = restPos.row(v[3]) - restPos.row(v[1]);
            double c01 = cotan(e0, e1);
            double c02 = cotan(e0, e2);
            double c03 = cotan(-e0, e3);
            double c04 = cotan(-e0, e4);
            Eigen::Vector4d K0(
                c03 + c04,
                c01 + c02,
                -c01 - c03,
                -c02 - c04);
            double A0 = e0.cross(e1).norm() / 2.0;
            double A1 = e0.cross(e2).norm() / 2.0;
            Eigen::Matrix4d Q = 3.0 / (A0 + A1) * weight * K0 * K0.transpose();
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        if (Q(j, k) != 0.0)
                        {
                            Qentries.emplace_back(3 * v[j] + l, 3 * v[k] + l, Q(j, k));
                        }
                    }
                }
            }
        }
        bendingM_.resize(3 * nverts, 3 * nverts);
        bendingM_.setFromTriplets(Qentries.begin(), Qentries.end());
    }
};



#endif