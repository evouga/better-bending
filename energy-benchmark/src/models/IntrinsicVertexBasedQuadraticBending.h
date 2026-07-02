#ifndef INTRINSICVERTEXBASEDQUADRATICBENDING_H
#define INTRINSICVERTEXBASEDQUADRATICBENDING_H

#include "Model.h"
#include <Eigen/Sparse>
#include <igl/intrinsic_delaunay_cotmatrix.h>
#include <igl/edge_lengths.h>
#include <igl/intrinsic_delaunay_triangulation.h>

class IntrinsicVertexBasedQuadraticBending : public Model
{
private:
    void buildBendingMatrix(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& restPos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        Eigen::SparseMatrix<double> &bendingM)
    {
        int nverts = restPos.rows();
        int nfaces = mesh.nFaces();
        int nedges = mesh.nEdges();

        Eigen::SparseMatrix<double> L;
        igl::intrinsic_delaunay_cotmatrix(restPos, mesh.faces(), L);
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
            if(mesh.edgeFace(i, 0) == -1 || mesh.edgeFace(i, 1) == -1)
            {
                int v0idx = mesh.edgeVertex(i, 0);
                int v1idx = mesh.edgeVertex(i, 1);
                bdry[v0idx] = true;
                bdry[v1idx] = true;
			}
        }

        Eigen::MatrixXd edgeLengths;
		igl::edge_lengths(restPos, mesh.faces(), edgeLengths);

        Eigen::MatrixXi newFaces;
		Eigen::MatrixXd newLengths;
        igl::intrinsic_delaunay_triangulation(edgeLengths, mesh.faces(), newLengths, newFaces);
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
            Eigen::Vector3d e1 = restPos.row(mesh.faceVertex(i, 1)) - restPos.row(mesh.faceVertex(i, 0));
			Eigen::Vector3d e2 = restPos.row(mesh.faceVertex(i, 2)) - restPos.row(mesh.faceVertex(i, 0));
			double area = 0.5 * e1.cross(e2).norm();            

            for (int j = 0; j < 3; j++)
            {
                int vidx = mesh.faceVertex(i, j);
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

        Eigen::SparseMatrix<double> biL = (bigL.transpose() + N.transpose()) * Minv * (bigL + N);

        std::vector<Eigen::Triplet<double> > bendingMcoeffs_;
        for (int k = 0; k < biL.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(biL, k); it; ++it)
            {
                bendingMcoeffs_.push_back({ (int)it.row(), (int)it.col(), it.value() });
            }
        }


        bendingM.resize(3 * nverts, 3 * nverts);
        bendingM.setFromTriplets(bendingMcoeffs_.begin(), bendingMcoeffs_.end());
    }

    double elasticEnergy(const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& restPos,
        const Eigen::MatrixXd& curPos,
        double thickness,
        double lameAlpha,
        double lameBeta)
    {
        int nverts = restPos.rows();
        Eigen::VectorXd displacement(3 * nverts);
        for (int i = 0; i < nverts; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                displacement[3 * i + j] = curPos(i, j);
            }
        }

        Eigen::SparseMatrix<double> bendingM;
        buildBendingMatrix(mesh, restPos, thickness, lameAlpha, lameBeta, bendingM);
        return 0.5 * displacement.transpose() * bendingM * displacement;
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
		return elasticEnergy(mesh, restPos, curPos, thickness, lameAlpha, lameBeta);
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
        return elasticEnergy(mesh, restPos, curPos, thickness, lameAlpha, lameBeta);
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
        return elasticEnergy(mesh, restPos, curPos, thickness, lameAlpha, lameBeta);
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
        return elasticEnergy(mesh, curPos, curPos, thickness, lameAlpha, lameBeta);
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
        return elasticEnergy(mesh, curPos, curPos, thickness, lameAlpha, lameBeta);
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
        return elasticEnergy(mesh, restPos, curPos, thickness, lameAlpha, lameBeta);
    }
};

#endif