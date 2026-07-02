#ifndef QBPL_H
#define QBPL_H

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <igl/intrinsic_delaunay_cotmatrix.h>

#include "Model.h"
#include "../include/MeshConnectivity.h"
#include "../include/MidedgeAngleTanFormulation.h"
#include "../include/ElasticShell.h"
#include "../include/RestState.h"
#include "../include/StVKMaterial.h"

class QBPL : public Model
{
public:
    QBPL(const Eigen::MatrixXd& restV, const Eigen::MatrixXi& restF,
        double h, double lameAlpha, double lameBeta)
        : mesh_(restF), restPos_(restV), mat_()
    {
        buildQBPLBendingMatrix(restV, h, lameAlpha, lameBeta);

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

    void buildQBPLBendingMatrix(
        const Eigen::MatrixXd& restPos,
        double thickness,
        double lameAlpha,
        double lameBeta)
    {
        int nverts = restPos.rows();
        int nfaces = mesh_.nFaces();
        int nedges = mesh_.nEdges();

        Eigen::SparseMatrix<double> L;
        igl::intrinsic_delaunay_cotmatrix(restPos, mesh_.faces(), L);
        std::vector<Eigen::Triplet<double> > bigLcoeffs;
        for (int k = 0; k < L.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
            {
                for (int j = 0; j < 3; j++)
                {
                    bigLcoeffs.push_back({ 3 * (int)it.row() + j, 3 * (int)it.col() + j, it.value() });
                }
            }
        }

        Eigen::SparseMatrix<double> bigL(3 * nverts, 3 * nverts);
        bigL.setFromTriplets(bigLcoeffs.begin(), bigLcoeffs.end());

        std::vector<bool> bdry(nverts);
        for (int i = 0; i < nedges; i++)
        {
            if (mesh_.edgeFace(i, 0) == -1 || mesh_.edgeFace(i, 1) == -1)
            {
                int v0idx = mesh_.edgeVertex(i, 0);
                int v1idx = mesh_.edgeVertex(i, 1);
                bdry[v0idx] = true;
                bdry[v1idx] = true;
            }
        }

        Eigen::MatrixXd edgeLengths;
        igl::edge_lengths(restPos, mesh_.faces(), edgeLengths);

        Eigen::MatrixXi newFaces;
        Eigen::MatrixXd newLengths;
        igl::intrinsic_delaunay_triangulation(edgeLengths, mesh_.faces(), newLengths, newFaces);
        int nnewfaces = newFaces.rows();

        std::vector<Eigen::Triplet<double> > Ncoeffs;
        for (int i = 0; i < nnewfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int v0idx = newFaces(i, j);
                int v1idx = newFaces(i, (j + 1) % 3);
                int v2idx = newFaces(i, (j + 2) % 3);
                if (bdry[v0idx] && bdry[v1idx])
                {
                    double l0 = newLengths(i, j);
                    double l1 = newLengths(i, (j + 1) % 3);
                    double l2 = newLengths(i, (j + 2) % 3);
                    double S = 0.25 * sqrt((l0 + l1 + l2) * (l0 + l1 - l2) * (l0 - l1 + l2) * (-l0 + l1 + l2));
                    double cot0 = (l1 * l1 + l2 * l2 - l0 * l0) / (4.0 * S);
                    double cot1 = (l0 * l0 + l2 * l2 - l1 * l1) / (4.0 * S);
                    for (int k = 0; k < 3; k++)
                    {
                        Ncoeffs.push_back({ 3 * v0idx + k, 3 * v0idx + k, 0.5 * cot1 });
                        Ncoeffs.push_back({ 3 * v0idx + k, 3 * v1idx + k, 0.5 * cot0 });
                        Ncoeffs.push_back({ 3 * v0idx + k, 3 * v2idx + k, 0.5 * (-cot0 - cot1) });
                        Ncoeffs.push_back({ 3 * v1idx + k, 3 * v0idx + k, 0.5 * cot1 });
                        Ncoeffs.push_back({ 3 * v1idx + k, 3 * v1idx + k, 0.5 * cot0 });
                        Ncoeffs.push_back({ 3 * v1idx + k, 3 * v2idx + k, 0.5 * (-cot0 - cot1) });
                    }
                }
            }
        }

        Eigen::SparseMatrix<double> N(3 * nverts, 3 * nverts);
        N.setFromTriplets(Ncoeffs.begin(), Ncoeffs.end());

        std::vector<double> Mcoeffs(nverts);
        std::vector<double> energycoeffs(nverts);
        std::vector<Eigen::Triplet<double> > Minvcoeffs;

        for (int i = 0; i < nfaces; i++)
        {
            double weight = thickness * thickness * thickness / 12.0 * (lameAlpha + 2.0 * lameBeta);
            Eigen::Vector3d e1 = restPos.row(mesh_.faceVertex(i, 1)) - restPos.row(mesh_.faceVertex(i, 0));
            Eigen::Vector3d e2 = restPos.row(mesh_.faceVertex(i, 2)) - restPos.row(mesh_.faceVertex(i, 0));
            double area = 0.5 * e1.cross(e2).norm();

            for (int j = 0; j < 3; j++)
            {
                int vidx = mesh_.faceVertex(i, j);
                Mcoeffs[vidx] += area / 3.0;
                energycoeffs[vidx] += weight * area / 3.0;
            }
        }
        for (int i = 0; i < nverts; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Minvcoeffs.push_back({ 3 * i + j, 3 * i + j, energycoeffs[i] / Mcoeffs[i] / Mcoeffs[i] });
            }
        }
        Eigen::SparseMatrix<double> Minv(3 * nverts, 3 * nverts);
        Minv.setFromTriplets(Minvcoeffs.begin(), Minvcoeffs.end());

        bendingM_ = (bigL.transpose() + N.transpose()) * Minv * (bigL + N);
    }
};



#endif