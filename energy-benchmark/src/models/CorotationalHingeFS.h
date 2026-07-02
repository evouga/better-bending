#ifndef COROTATIONALHINGEFS_H
#define COROTATIONALHINGEFS_H

#include "Model.h"
#include <Eigen/Sparse>

class CorotationalHingeFS : public Model
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
        int nfaces = mesh.nFaces();

		double weight = thickness * thickness * thickness / 12.0 * (lameAlpha + 2.0 * lameBeta);
        double result = 0;

        for (int i = 0; i < nfaces; i++)
        {
            const Eigen::MatrixXd* refPos = restPos ? restPos : &curPos;

            Eigen::Vector3d v0bar = refPos->row(mesh.faceVertex(i, 0));
            Eigen::Vector3d v1bar = refPos->row(mesh.faceVertex(i, 1));
            Eigen::Vector3d v2bar = refPos->row(mesh.faceVertex(i, 2));
            Eigen::Vector3d nbar = (v1bar - v0bar).cross(v2bar - v0bar);
			Eigen::Vector3d nbarhat = nbar.normalized();
            double Abar = 0.5 * nbar.norm();
            Eigen::Vector3d ebar[3];
            ebar[0] = v2bar - v1bar;
            ebar[1] = v0bar - v2bar;
            ebar[2] = v1bar - v0bar;
            
            Eigen::Matrix3d RTDR;
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    RTDR(j, k) = ebar[j].dot(ebar[k]) * ebar[j].dot(ebar[k]) / ebar[j].squaredNorm() / ebar[k].squaredNorm();
                }
            }

			Eigen::Vector3d v0 = curPos.row(mesh.faceVertex(i, 0));
			Eigen::Vector3d v1 = curPos.row(mesh.faceVertex(i, 1));
            Eigen::Vector3d v2 = curPos.row(mesh.faceVertex(i, 2));
			Eigen::Vector3d nhat = (v1 - v0).cross(v2 - v0).normalized();

            Eigen::Vector3d LpX;
            Eigen::Vector3d Lpx;
            for (int j = 0; j < 3; j++)
            {
                int oppv = mesh.vertexOppositeFaceEdge(i, j);
                if (oppv == -1)
                {
                    LpX[j] = 0;
                    Lpx[j] = 0;
                }
                else
                {
                    Eigen::Vector3d vOpp = refPos->row(oppv);
					Eigen::Vector3d vp1 = refPos->row(mesh.faceVertex(i, (j + 1) % 3));
                    Eigen::Vector3d vp2 = refPos->row(mesh.faceVertex(i, (j + 2) % 3));
					double Aopp = 0.5 * (vp1 - vOpp).cross(vp2 - vOpp).norm();

					Eigen::Vector3d vOppCur = curPos.row(oppv);
                    double X = (vOppCur - v0).dot(nhat);
                    double x = 0;
                    if (restPos)
                    {
						x = (vOpp - v0bar).dot(nbarhat);
                    }

                    LpX[j] = X * ebar[j].squaredNorm() / (Abar + Aopp) / 2.0 / Aopp;
					Lpx[j] = x * ebar[j].squaredNorm() / (Abar + Aopp) / 2.0 / Aopp;
                }
            }

			result += 0.5 * Abar * weight * (LpX - Lpx).transpose() * RTDR * (LpX - Lpx);
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