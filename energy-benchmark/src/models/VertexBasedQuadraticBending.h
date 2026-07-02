#ifndef VERTEXBASEDQUADRATICBENDING_H
#define VERTEXBASEDQUADRATICBENDING_H

#include "Model.h"
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>

class VertexBasedQuadraticBending : public Model
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
        igl::cotmatrix(restPos, mesh.faces(), L);
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
        std::vector<Eigen::Triplet<double> > Ncoeffs;
        for (int i = 0; i < nedges; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                if (mesh.edgeFace(i, 1 - side) == -1)
                {
                    int oppvidx = mesh.edgeOppositeVertex(i, side);
                    int v0idx = mesh.edgeVertex(i, 0);
                    int v1idx = mesh.edgeVertex(i, 1);
                    bdry[v0idx] = true;
                    bdry[v1idx] = true;
                    Eigen::Vector3d oppv = restPos.row(oppvidx).transpose();
                    Eigen::Vector3d v0 = restPos.row(v0idx).transpose();
                    Eigen::Vector3d v1 = restPos.row(v1idx).transpose();
                    double cot0 = (oppv - v0).dot(v1 - v0) / (oppv - v0).cross(v1 - v0).norm();
                    double cot1 = (oppv - v1).dot(v0 - v1) / (oppv - v1).cross(v0 - v1).norm();
                    for (int j = 0; j < 3; j++)
                    {
                        Ncoeffs.push_back({ 3 * v0idx + j, 3 * v0idx + j, 0.5 * cot1 });
                        Ncoeffs.push_back({ 3 * v0idx + j, 3 * v1idx + j, 0.5 * cot0 });
                        Ncoeffs.push_back({ 3 * v0idx + j, 3 * oppvidx + j, 0.5 * (-cot0 - cot1) });
                        Ncoeffs.push_back({ 3 * v1idx + j, 3 * v0idx + j, 0.5 * cot1 });
                        Ncoeffs.push_back({ 3 * v1idx + j, 3 * v1idx + j, 0.5 * cot0 });
                        Ncoeffs.push_back({ 3 * v1idx + j, 3 * oppvidx + j, 0.5 * (-cot0 - cot1) });
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