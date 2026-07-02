#ifndef DISCRETESHELLSTAN_H
#define DISCRETESHELLSTAN_H

#include "Model.h"
#include <Eigen/Sparse>

class DiscreteShellsTan : public Model
{
private:
    double elasticEnergy(const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd *restPos,
        const Eigen::MatrixXd& curPos,
        double thickness,
        double lameAlpha,
        double lameBeta)
    {
        int nverts = curPos.rows();
		int nedges = mesh.nEdges();
		double weight = thickness * thickness * thickness / 12.0 * (lameAlpha + 2.0 * lameBeta);
        double result = 0;
        for (int i = 0; i < nedges; i++)
        {
            if(mesh.edgeFace(i, 0) == -1 || mesh.edgeFace(i, 1) == -1)
				continue; // skip boundary edges
            int v[4];
			v[0] = mesh.edgeVertex(i, 0);
			v[1] = mesh.edgeVertex(i, 1);
			v[2] = mesh.edgeOppositeVertex(i, 0);
            v[3] = mesh.edgeOppositeVertex(i, 1);
            double thetabar = 0;
            double coeff = 0;

            Eigen::Vector3d e0 = curPos.row(v[1]) - curPos.row(v[0]);
            Eigen::Vector3d e1 = curPos.row(v[2]) - curPos.row(v[0]);
            Eigen::Vector3d e2 = curPos.row(v[3]) - curPos.row(v[0]);
            Eigen::Vector3d n0 = e1.cross(e0).normalized();
            Eigen::Vector3d n1 = e0.cross(e2).normalized();
            double theta = std::acos(std::min(1.0, n0.dot(n1)));
            if (restPos)
            {
				Eigen::Vector3d e0bar = restPos->row(v[1]) - restPos->row(v[0]);
                Eigen::Vector3d e1bar = restPos->row(v[2]) - restPos->row(v[0]);
				Eigen::Vector3d e2bar = restPos->row(v[3]) - restPos->row(v[0]);
                Eigen::Vector3d n0bar = e1bar.cross(e0bar).normalized();
                Eigen::Vector3d n1bar = e0bar.cross(e2bar).normalized();
				thetabar = std::acos(std::min(1.0, n0bar.dot(n1bar)));
				double A1bar = e0bar.cross(e1bar).norm() / 2.0;
				double A2bar = e0bar.cross(e2bar).norm() / 2.0;
                coeff = 3.0 * e0bar.squaredNorm() / (A1bar + A2bar) * weight;
            }
            else
			{
				double A1 = e0.cross(e1).norm() / 2.0;
				double A2 = e0.cross(e2).norm() / 2.0;
				coeff = 3.0 * e0.squaredNorm() / (A1 + A2) * weight;
            }
            result += 0.5 * coeff * (2.0 * std::tan(theta / 2.0) - 2.0 * std::tan(thetabar / 2.0)) * (2.0 * std::tan(theta / 2.0) - 2.0 * std::tan(thetabar / 2.0));
        }
        return result;
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
        const Eigen::MatrixXd &curEdgePos,
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