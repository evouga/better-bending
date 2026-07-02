#ifndef SIXVERTEXMIDEDGEAVERAGE_H
#define SIXVERTEXMIDEDGEAVERAGE_H

#include "Model.h"
#include "../include/MidedgeAverageFormulation.h"
#include "../include/RestState.h"
#include "../include/ElasticShell.h"
#include "../include/StVKMaterial.h"

class SixVertexMidedgeAverage : public Model
{
private:
    double elasticEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd* restPos,
        const Eigen::MatrixXd& curPos,
        double thickness,
        double lameAlpha,
        double lameBeta)
    {
        Eigen::VectorXd edgeDOFs;
		LibShell::MidedgeAverageFormulation::initializeExtraDOFs(edgeDOFs, mesh, restPos ? *restPos : curPos);            

        // initialize the rest geometry of the shell
        LibShell::MonolayerRestState restState;

        // set uniform thicknesses
        restState.thicknesses.resize(mesh.nFaces(), thickness);
        restState.lameAlpha.resize(mesh.nFaces(), lameAlpha);
        restState.lameBeta.resize(mesh.nFaces(), lameBeta);

        // initialize first and second fundamental forms to those of input mesh
        LibShell::ElasticShell<LibShell::MidedgeAverageFormulation>::
            firstFundamentalForms(mesh, restPos ? *restPos : curPos, restState.abars);
        if (restPos)
        {
            LibShell::ElasticShell<LibShell::MidedgeAverageFormulation>::
                secondFundamentalForms(mesh, *restPos, edgeDOFs, restState.bbars);
        }
        else
        {
            restState.bbars.resize(mesh.nFaces());
            for(int i = 0; i < mesh.nFaces(); i++)
            {
                restState.bbars[i].setZero();
			}
        }

        LibShell::StVKMaterial<LibShell::MidedgeAverageFormulation> mat;

        return LibShell::ElasticShell<LibShell::MidedgeAverageFormulation>::elasticEnergy(mesh, curPos, edgeDOFs, mat, restState, NULL, NULL);
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
