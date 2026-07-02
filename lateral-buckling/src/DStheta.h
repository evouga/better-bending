#ifndef DSTHETA_H
#define DSTHETA_H

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>

#include <igl/intrinsic_delaunay_cotmatrix.h>

#include "Model.h"
#include "../include/MeshConnectivity.h"
#include "../include/MidedgeAngleTanFormulation.h"
#include "../include/ElasticShell.h"
#include "../include/RestState.h"
#include "../include/StVKMaterial.h"
#include "../src/GeometryDerivatives.h"

class DSTheta : public Model
{
public:
    DSTheta(const Eigen::MatrixXd& restV, const Eigen::MatrixXi& restF,
        double h, double lameAlpha, double lameBeta)
        : mesh_(restF), restPos_(restV), mat_()
    {
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

        matCoeff_ = h * h * h / 12.0 * (lameAlpha + 2.0 * lameBeta);
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

        int nedges = mesh_.nEdges();
        for (int i = 0; i < nedges; i++)
        {
            if (mesh_.edgeFace(i, 0) == -1 || mesh_.edgeFace(i, 1) == -1)
                continue; // skip boundary edges
            
			int f1 = mesh_.edgeFace(i, 0);
			int f2 = mesh_.edgeFace(i, 1);

			int v1 = mesh_.edgeVertex(i, 0);
            int v2 = mesh_.edgeVertex(i, 1);

			Eigen::Vector3d e = curPos.row(v2) - curPos.row(v1);
			Eigen::Vector3d ebar = restPos_.row(v2) - restPos_.row(v1);

            Eigen::Matrix<double, 3, 9> n1deriv;
            std::vector<Eigen::Matrix<double, 9, 9> > n1hess;
            Eigen::Vector3d n1 = LibShell::faceNormal(mesh_, curPos, f1, 0, &n1deriv, &n1hess);

            Eigen::Matrix<double, 3, 9> n2deriv;
            std::vector<Eigen::Matrix<double, 9, 9> > n2hess;
            Eigen::Vector3d n2 = LibShell::faceNormal(mesh_, curPos, f2, 0, &n2deriv, &n2hess);

            Eigen::Vector3d n1bar = LibShell::faceNormal(mesh_, restPos_, f1, 0, nullptr, nullptr);
			Eigen::Vector3d n2bar = LibShell::faceNormal(mesh_, restPos_, f2, 0, nullptr, nullptr);
			
            double thetabar = LibShell::angle(n1bar, n2bar, ebar, nullptr, nullptr);

			Eigen::Matrix<double, 1, 9> thetaderiv;
            Eigen::Matrix<double, 9, 9> thetahess;
            double theta = LibShell::angle(n1, n2, e, &thetaderiv, &thetahess);

            double A1bar = n1bar.norm() / 2.0;
            double A2bar = n2bar.norm() / 2.0;
            double coeff = 3.0 * ebar.squaredNorm() / (A1bar + A2bar) * matCoeff_;

            energy += 0.5 * coeff * (theta - thetabar) * (theta - thetabar);

            if (deriv)
            {
				Eigen::Matrix<double, 1, 9> dElocal = coeff * (theta - thetabar) * thetaderiv;
                Eigen::Matrix<double, 1, 9> dEdf1 = dElocal.segment<3>(0) * n1deriv;
                Eigen::Matrix<double, 1, 9> dEdf2 = dElocal.segment<3>(3) * n2deriv;
                for (int j = 0; j < 3; j ++)
                {
					D.segment(3 * mesh_.faceVertex(f1, j), 3) += dEdf1.segment<3>(3 * j);
                    D.segment(3 * mesh_.faceVertex(f2, j), 3) += dEdf2.segment<3>(3 * j);
                }
            }

            if (hess)
            {
				Eigen::Matrix<double, 9, 9> d2Elocal = coeff * (thetaderiv.transpose() * thetaderiv + (theta - thetabar) * thetahess);

				Eigen::Matrix<double, 9, 9> d2Edf1df1 = n1deriv.transpose() * d2Elocal.block(0, 0, 3, 3) * n1deriv;
                Eigen::Matrix<double, 9, 9> d2Edf1df2 = n1deriv.transpose() * d2Elocal.block(0, 3, 3, 3) * n2deriv;
                Eigen::Matrix<double, 9, 3> d2Edf1de = n1deriv.transpose() * d2Elocal.block(0, 6, 3, 3);
                Eigen::Matrix<double, 9, 9> d2Edf2df1 = n2deriv.transpose() * d2Elocal.block(3, 0, 3, 3) * n1deriv;
                Eigen::Matrix<double, 9, 9> d2Edf2df2 = n2deriv.transpose() * d2Elocal.block(3, 3, 3, 3) * n2deriv;
                Eigen::Matrix<double, 9, 3> d2Edf2de = n2deriv.transpose() * d2Elocal.block(3, 6, 3, 3);
                Eigen::Matrix<double, 3, 9> d2Ededf1 = d2Elocal.block(6, 0, 3, 3) * n1deriv;
                Eigen::Matrix<double, 3, 9> d2Ededf2 = d2Elocal.block(6, 3, 3, 3) * n2deriv;
                Eigen::Matrix<double, 3, 3> d2Edede = d2Elocal.block(6, 6, 3, 3);
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            for (int m = 0; m < 3; m++)
                            {
                                Hcoeffs.push_back({ 3 * mesh_.faceVertex(f1, j) + l, 3 * mesh_.faceVertex(f1, k) + m, d2Edf1df1(3 * j + l, 3 * k + m) });
                                Hcoeffs.push_back({ 3 * mesh_.faceVertex(f1, j) + l, 3 * mesh_.faceVertex(f2, k) + m, d2Edf1df2(3 * j + l, 3 * k + m) });
                                Hcoeffs.push_back({ 3 * mesh_.faceVertex(f2, j) + l, 3 * mesh_.faceVertex(f1, k) + m, d2Edf2df1(3 * j + l, 3 * k + m) });
                                Hcoeffs.push_back({ 3 * mesh_.faceVertex(f2, j) + l, 3 * mesh_.faceVertex(f2, k) + m, d2Edf2df2(3 * j + l, 3 * k + m) });
                            }
							Hcoeffs.push_back({ 3 * mesh_.faceVertex(f1, j) + l, 3 * v1 + k, -d2Edf1de(3 * j + l, k) });
                            Hcoeffs.push_back({ 3 * mesh_.faceVertex(f1, j) + l, 3 * v2 + k, d2Edf1de(3 * j + l, k) });
                            Hcoeffs.push_back({ 3 * mesh_.faceVertex(f2, j) + l, 3 * v1 + k, -d2Edf2de(3 * j + l, k) });
                            Hcoeffs.push_back({ 3 * mesh_.faceVertex(f2, j) + l, 3 * v2 + k, d2Edf2de(3 * j + l, k) });
                        }						
						Hcoeffs.push_back({ 3 * v1 + j, 3 * v1 + k, d2Edede(j, k) });
                        Hcoeffs.push_back({ 3 * v1 + j, 3 * v2 + k, -d2Edede(j, k) });
                        Hcoeffs.push_back({ 3 * v2 + j, 3 * v1 + k, -d2Edede(j, k) });
						Hcoeffs.push_back({ 3 * v2 + j, 3 * v2 + k, d2Edede(j, k) });
                    }
                }
            }
        }

        if (deriv)
        {
            *deriv = D.segment(0, nposdofs);
        }

        if (hess)
        {
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
    LibShell::StVKMaterial<LibShell::MidedgeAngleTanFormulation> mat_;
    LibShell::MonolayerRestState restState_;
    double matCoeff_;
};



#endif
