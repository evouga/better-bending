#ifndef COROTATIONALHINGEES_H
#define COROTATIONALHINGEES_H

#include "Model.h"
#include <Eigen/Sparse>

class CorotationalHingeES : public Model
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
            
            double Akb = 0;

            Eigen::Vector3d e0 = curPos.row(v[1]) - curPos.row(v[0]);
            Eigen::Vector3d e1 = curPos.row(v[2]) - curPos.row(v[0]);
            Eigen::Vector3d e2 = curPos.row(v[3]) - curPos.row(v[0]);
            Eigen::Vector3d n0 = e1.cross(e0).normalized();
            Eigen::Vector3d n1 = e0.cross(e2).normalized();
			Eigen::Vector3d n = (n0 + n1).normalized();

            double km0 = 0;
            double km = 0;

            if (restPos)
            {
				Eigen::Vector3d e0bar = restPos->row(v[1]) - restPos->row(v[0]);
                Eigen::Vector3d e1bar = restPos->row(v[2]) - restPos->row(v[0]);
				Eigen::Vector3d e2bar = restPos->row(v[3]) - restPos->row(v[0]);
                Eigen::Vector3d n0bar = e1bar.cross(e0bar).normalized();
                Eigen::Vector3d n1bar = e0bar.cross(e2bar).normalized();
				
                double A1bar = e0bar.cross(e1bar).norm() / 2.0;
				double A2bar = e0bar.cross(e2bar).norm() / 2.0;
                
				Eigen::Vector3d nbar = (n0bar + n1bar).normalized();

				double Lmcoeff = e0bar.norm() / (A1bar + A2bar);
                km0 = Lmcoeff * e0bar.norm() * (nbar.dot(e1bar) / 2.0 / A1bar + nbar.dot(e2bar) / 2.0 / A2bar);
                km = Lmcoeff * e0bar.norm() * (n.dot(e1) / 2.0 / A1bar + n.dot(e2) / 2.0 / A2bar);

                Akb = (A1bar + A2bar) * weight;
            }
            else
			{
				double A1 = e0.cross(e1).norm() / 2.0;
				double A2 = e0.cross(e2).norm() / 2.0;
				
                double Lmcoeff = e0.norm() / (A1 + A2);

                km = Lmcoeff * e0.norm() * (n.dot(e1) / 2.0 / A1 + n.dot(e2) / 2.0 / A2);

                Akb = (A1 + A2) * weight;
            }
            result += 0.5 * Akb * (km - km0) * (km - km0);
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