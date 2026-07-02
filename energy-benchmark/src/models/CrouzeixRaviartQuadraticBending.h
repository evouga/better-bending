#ifndef CROUZEIXRAVIARTQUADRATICBENDING_H
#define CROUZEIXRAVIARTQUADRATICBENDING_H

#include "Model.h"
#include <Eigen/Sparse>

static double cotan(Eigen::Vector3d e0, Eigen::Vector3d e1)
{
    double dot = e0.dot(e1);
    double crossNorm = e0.cross(e1).norm();
    if (crossNorm == 0.0)
        return 0.0;
    return dot / crossNorm;
}

class CrouzeixRaviartQuadraticBending : public Model
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

        double weight = thickness * thickness * thickness / 12.0 * (lameAlpha + 2.0 * lameBeta);
        std::vector<Eigen::Triplet<double>> Qentries;
        for (int i = 0; i < nedges; i++)
        {
            if(mesh.edgeOppositeVertex(i, 0) == -1 || mesh.edgeOppositeVertex(i, 1) == -1)
				continue; // skip boundary edges    
            int v[4];
			v[0] = mesh.edgeVertex(i, 0);
			v[1] = mesh.edgeVertex(i, 1);
			v[2] = mesh.edgeOppositeVertex(i, 0);
			v[3] = mesh.edgeOppositeVertex(i, 1);
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
		bendingM.resize(3 * nverts, 3 * nverts);
        bendingM.setFromTriplets(Qentries.begin(), Qentries.end());        
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
        const Eigen::MatrixXd &restEdgePos,
        const Eigen::MatrixXd& curPos,
        const Eigen::MatrixXd &curEdgePos,
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
        const Eigen::MatrixXd &curEdgePos,
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
        const Eigen::MatrixXd &curEdgePos,
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
        const Eigen::MatrixXd &curEdgePos,
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