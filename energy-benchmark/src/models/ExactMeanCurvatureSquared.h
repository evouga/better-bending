#ifndef EXACTMEANCURVATURESQUARED_H
#define EXACTMEANCURVATURESQUARED_H

#include "Model.h"

// Analytic "mean curvature style" elastic energy

class ExactMeanCurvatureSquared : public Model
{
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
        // ground truth energy
        // W = PI * r
        // r(x,y) = (r cos[x/r], r sin[x/r], y)^T
        // dr(x,y) = ((-sin[x/r], 0),
        //            ( cos[x/r], 0),
        //            ( 0, 1 ))
        Eigen::Matrix2d abar;
        abar.setIdentity();

        // n = (-sin[x/r], cos[x/r], 0) x (0, 0, 1) = ( cos[x/r], sin[x/r], 0 )
        // dn = ((-sin[x/r]/r, 0),
        //       ( cos[x/r]/r, 0),
        //       ( 0, 0 ))
        // b = dr^T dn = ((1/r, 0), (0, 0))
        Eigen::Matrix2d b;
        b << 1.0 / radius, 0, 0, 0;

        Eigen::Matrix2d M = abar.inverse() * b;
        double meancurvnorm = (lameAlpha / 2.0 + lameBeta) * M.trace() * M.trace();
        double coeff = thickness * thickness * thickness / 12.0;
        constexpr double PI = 3.1415926535898;
        double area = PI * radius * height;

        return meancurvnorm * coeff * area;
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
        // ground truth energy
        // W = PI * r
        // r(x,y) = (r cos[x/r], r sin[x/r], y)^T
        // dr(x,y) = ((-sin[x/r], 0),
        //            ( cos[x/r], 0),
        //            ( 0, 1 ))
        Eigen::Matrix2d abar;
        abar.setIdentity();

        // n = (-sin[x/r], cos[x/r], 0) x (0, 0, 1) = ( cos[x/r], sin[x/r], 0 )
        // dn = ((-sin[x/r]/r, 0),
        //       ( cos[x/r]/r, 0),
        //       ( 0, 0 ))
        // b = dr^T dn = ((1/r, 0), (0, 0))
        Eigen::Matrix2d b;
        b << 1.0 / radius, 0, 0, -1.0 / radius;

        Eigen::Matrix2d M = abar.inverse() * b;
        double meancurvnorm = (lameAlpha / 2.0 + lameBeta) * M.trace() * M.trace();
        double coeff = thickness * thickness * thickness / 12.0;
        constexpr double PI = 3.1415926535898;
        double area = PI * radius * height;

        return meancurvnorm * coeff * area;
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
		// crease has infinite curvature, so infinite energy
		return std::numeric_limits<double>::infinity();
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
        Eigen::Matrix2d abar;
        abar.setIdentity();

        Eigen::Matrix2d b;
        b << 1.0 / radius, 0, 0, 1.0 / radius;

        Eigen::Matrix2d M = abar.inverse() * b;
        double meancurvnorm = (lameAlpha / 2.0 + lameBeta) * M.trace() * M.trace();
        double coeff = thickness * thickness * thickness / 12.0;
        constexpr double PI = 3.1415926535898;
        double area = 4.0 * PI * radius * radius;

        return meancurvnorm * coeff * area;
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
        // ground truth energy
        // r(x,y) = (R1 + R2 cos[y]) cos[x], (R1 + R2 cos[y]) sin[x], R2 sin[y])^T
        // dr(x,y) = (( -(R1 + R2 cos[y]) sin[x], -R2 sin[y] cos[x]),
        //            ( (R1 + R2 cos[y]) cos[x], -R2 sin[y] sin[x]),
        //            ( 0, R2 cos[y]))
        // a(x,y) = (( (R1 + R2 cos[y])^2, 0)
        //           ( 0, R2^2) )
        // da = R2 (R1 + R2 cos[y])
        // b = (( (R1 + R2 cos[y]) cos[y], 0)
        //      ( 0, R2 ))


        double coeff = thickness * thickness * thickness / 12.0;
        constexpr double PI = 3.1415926535898;
        double sliceW = PI * outerRadius * outerRadius * (lameAlpha + 2.0 * lameBeta) / 2.0 / innerRadius / std::sqrt(outerRadius * outerRadius - innerRadius * innerRadius);
        return coeff * sliceW * 2.0 * PI;
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
        double coeff = thickness * thickness * thickness / 12.0;
        constexpr double PI = 3.1415926535898;
        double num = helixHeight * helixHeight * (lameAlpha + 2.0 * lameBeta) * std::log(2.0);
        double denom = 2.0 * helixRadius * std::sqrt(helixHeight*helixHeight + 4.0 * PI * PI * helixRadius * helixRadius);
        return coeff * num / denom;
    }
};

#endif
