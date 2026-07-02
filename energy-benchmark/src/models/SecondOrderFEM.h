#ifndef SECONDORDERFEM_H
#define SECONDORDERFEM_H

#include "Model.h"
#include <Eigen/Sparse>
#include "../Derivatives.h"
#include "../src/GeometryDerivatives.h"
#include "../Visualization.h"

struct QuadraticElement
{
    static Eigen::Vector<double, 6> N(double s, double t)
    {
        Eigen::Vector<double, 6> result;
        result << 2.0 * (1.0 - s - t) * (0.5 - s - t),
            s* (2.0 * s - 1.0),
            t* (2.0 * t - 1.0),
            4.0 * s * t,
            4.0 * t * (1.0 - s - t),
            4.0 * s * (1.0 - s - t);
        return result;
    }

    static Eigen::Matrix<double, 6, 2> dN(double s, double t)
    {
        Eigen::Matrix<double, 6, 2> result;
        result << -3.0 + 4.0 * s + 4.0 * t, -3.0 + 4.0 * s + 4.0 * t,
            4.0 * s - 1.0, 0.0,
            0, 4.0 * t - 1.0,
            4.0 * t, 4.0 * s,
            -4.0 * t, 4.0 * (1.0 - s - 2.0 * t),
            4.0 * (1.0 - 2.0 * s - t), -4.0 * s;
        return result;
    }

    static Eigen::Matrix<double, 6, 3> d2N(double s, double t)
    {
        Eigen::Matrix<double, 6, 3> result;
        result << 4.0, 4.0, 4.0,
            4.0, 0.0, 0.0,
            0.0, 0.0, 4.0,
            0.0, 4.0, 0.0,
            0.0, -4.0, -8.0,
            -8.0, -4.0, 0.0;
        return result;
    }
};

class Quadrature
{
public:
    Quadrature() : quadpoints(0) {}
    int numQuadraturePoint()
    {
        return quadpoints;
    }

    double s(int quadPt)
    {
		return svals[quadPt];
    }
    
    double t(int quadPt)
    {
        return tvals[quadPt];
    }

    double weight(int quadPt)
    {
        return wvals[quadPt];
    }

protected:
    int quadpoints;
    std::vector<double> svals;
    std::vector<double> tvals;
    std::vector<double> wvals;
};

struct SecondOrderEdgeQuadrature : public Quadrature
{
    SecondOrderEdgeQuadrature()
    {
        quadpoints = 2;
        svals = { 0.5 - 0.5 * std::sqrt(1.0 / 3.0), 0.5 + 0.5 * std::sqrt(1.0 / 3.0) };
        tvals = { 0.0, 0.0 };
        wvals = { 1.0 / 2.0, 1.0 / 2.0 };
    }
};

struct SecondOrderTriangleQuadrature : public Quadrature
{
    SecondOrderTriangleQuadrature()
    {
        quadpoints = 3;
        svals = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
        tvals = { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 };
        wvals = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
    }
};

struct SecondOrderEdgeBasedTriangleQuadrature : public Quadrature
{
    SecondOrderEdgeBasedTriangleQuadrature()
    {
        quadpoints = 7;
        svals = { 0.5 - 0.5 * std::sqrt(1.0 / 3.0), 0.5 + 0.5 * std::sqrt(1.0 / 3.0), 0.0, 0.0, 0.5 - 0.5 * std::sqrt(1.0 / 3.0), 0.5 + 0.5 * std::sqrt(1.0 / 3.0), 1.0 / 3.0 };
        tvals = { 0.0, 0.0, 0.5 - 0.5 * std::sqrt(1.0 / 3.0), 0.5 + 0.5 * std::sqrt(1.0 / 3.0), 0.5 + 0.5 * std::sqrt(1.0 / 3.0), 0.5 - 0.5 * std::sqrt(1.0 / 3.0), 1.0 / 3.0 };
        wvals = { 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 2.0 };
    }
};

struct FourthOrderEdgeBasedTriangleQuadrature : public Quadrature
{
    FourthOrderEdgeBasedTriangleQuadrature()
    {
        quadpoints = 13;
        double a = 0.5 - 0.5 * std::sqrt(1.0 / 3.0);
        svals = { a, 1.0 - a, 0.0, 0.0, a, 1.0 - a, 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 2.0, 0.0, 1.0 / 2.0, 1.0/3.0 };
        tvals = { 0.0, 0.0, a, 1.0 - a, 1.0 - a, a, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0, 0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0/3.0 };
        wvals = { 1.0 / 20.0, 1.0 / 20.0, 1.0 / 20.0, 1.0 / 20.0, 1.0 / 20.0, 1.0 / 20.0, 1.0 / 10.0, 1.0 / 10.0, 1.0 / 10.0, 1.0 / 30.0, 1.0 / 30.0, 1.0 / 30.0, 3.0 / 10.0 };
    }
};

class SecondOrderFEM : public Model
{
private:
    void visualize(const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& paramDomain,
        const Eigen::MatrixXd& vertPos,
        const Eigen::MatrixXd& edgePos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        Eigen::VectorXd& vertEnergyDensity,
        Eigen::MatrixXd& framePts,
        Eigen::MatrixXd& frameUs,
		Eigen::MatrixXd& frameVs,
        Eigen::MatrixXd& frameNs)
    {
		int nfaces = mesh.nFaces();
        constexpr int ptsPerEdge = 30;
        int ptsPerTriangle = ptsPerEdge * (ptsPerEdge + 1) / 2;
		int facesPerTriangle = (ptsPerEdge - 1) * (ptsPerEdge - 1);
		V.resize(nfaces * ptsPerTriangle, 3);
		F.resize(nfaces * facesPerTriangle, 3);
		vertEnergyDensity.resize(nfaces * ptsPerTriangle);
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < ptsPerEdge; j++)
            {
                for (int k = 0; k <= j; k++)
                {
					double s = double(j - k) / double(ptsPerEdge - 1);
					double t = double(k) / double(ptsPerEdge - 1);
					Eigen::Vector<double, 6> N = QuadraticElement::N(s, t);
                    Eigen::Matrix<double, 6, 2> dN = QuadraticElement::dN(s, t);
                    Eigen::Matrix<double, 6, 3> d2N = QuadraticElement::d2N(s, t);

					int idx = i * ptsPerTriangle + (j * (j + 1)) / 2 + k;
                    Eigen::Matrix<double, 18, 1> v;
                    for (int l = 0; l < 3; l++)
                    {
                        int vertIndex = mesh.faceVertex(i, l);
                        v.segment<3>(3 * l) = vertPos.row(vertIndex).transpose();
                        int edgeIndex = mesh.faceEdge(i, l);
                        v.segment<3>(9 + 3 * l) = edgePos.row(edgeIndex).transpose();                        
					}
                    Eigen::Vector3d pos(0, 0, 0);
                    for(int l=0; l< 6; l++)
						pos += N(l) * v.segment<3>(3 * l);
					V.row(idx) = pos.transpose();

                    Eigen::Matrix<double, 6, 18> dxdr;
                    dxdr.setZero();
                    for (int m = 0; m < 3; m++)
                    {
                        for (int l = 0; l < 6; l++)
                        {
                            dxdr(m, 3 * l + m) += dN(l, 0);
                            dxdr(3 + m, 3 * l + m) += dN(l, 1);
                        }
                    }

                    Eigen::Vector4d abar;
                    Eigen::Vector3d parame0 = paramDomain.row(mesh.faceVertex(i, 1)) - paramDomain.row(mesh.faceVertex(i, 0));
                    Eigen::Vector3d parame1 = paramDomain.row(mesh.faceVertex(i, 2)) - paramDomain.row(mesh.faceVertex(i, 0));

                    abar << parame0.dot(parame0), parame0.dot(parame1), parame0.dot(parame1), parame1.dot(parame1);

                    double dA = 0.5 * std::sqrt(abar[0] * abar[3] - abar[1] * abar[2]);

                    double coeff = thickness * thickness * thickness / 12.0;
                    Eigen::Vector4d bbar;
                    bbar.setZero();

                    std::vector<Eigen::Matrix<double, 3, 18> > dxd2r(4);
                    for (int k = 0; k < 4; k++)
                        dxd2r[k].setZero();
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 6; l++)
                        {
                            dxd2r[0](k, 3 * l + k) = d2N(l, 0);
                            dxd2r[1](k, 3 * l + k) = d2N(l, 1);
                            dxd2r[2](k, 3 * l + k) = d2N(l, 1);
                            dxd2r[3](k, 3 * l + k) = d2N(l, 2);
                        }
                    }

                    Eigen::Matrix<double, 6, 1> dr = dxdr * v;
                    std::vector<Eigen::Vector3d> d2r(4);
                    for (int k = 0; k < 4; k++)
                    {
                        d2r[k] = dxd2r[k] * v;
                    }

                    double meanCuvatureEnergy = MeanCurvatureSquaredEnergy<18>(dr, d2r, dxdr, dxd2r, abar, bbar, lameAlpha, lameBeta, NULL, NULL);
					vertEnergyDensity(idx) = coeff * meanCuvatureEnergy * dA;
                }
            }
            int cur = 0;
            for (int j = 0; j < ptsPerEdge - 1; j++)
            {
                for (int k = 0; k <= j; k++)
                {
                    int idx = i * ptsPerTriangle + (j * (j + 1)) / 2 + k;
                    int nextrow = i * ptsPerTriangle + ((j + 1) * (j + 2)) / 2 + k;
					int nextdiag = i * ptsPerTriangle + ((j + 1) * (j + 2)) / 2 + (k + 1);
					F(i * facesPerTriangle + cur, 0) = idx;
					F(i * facesPerTriangle + cur, 1) = nextrow;
					F(i * facesPerTriangle + cur, 2) = nextdiag;
					cur++;
                    if (k < j)
                    {
                        int nextcol = i * ptsPerTriangle + (j * (j + 1)) / 2 + k + 1;
                        F(i * facesPerTriangle + cur, 0) = idx;
                        F(i * facesPerTriangle + cur, 1) = nextdiag;
						F(i * facesPerTriangle + cur, 2) = nextcol;
                        cur++;
                    }
                }
            }
        }

        int interiorEdges = 0;
		int nedges = mesh.nEdges();
        for(int i=0; i<nedges; i++)
        {
            if(mesh.edgeFace(i, 0) != -1 && mesh.edgeFace(i, 1) != -1)
                interiorEdges++;
		}

		SecondOrderEdgeQuadrature edgeQuad;
		int quadpts = edgeQuad.numQuadraturePoint();
		framePts.resize(2 * interiorEdges * quadpts, 3);
		frameUs.resize(2 * interiorEdges * quadpts, 3);
		frameVs.resize(2 * interiorEdges * quadpts, 3);
		frameNs.resize(2 * interiorEdges * quadpts, 3);
        int idx = 0;
        for (int i = 0; i < nfaces; i++)
        {
            for (int k = 0; k < 3; k++)
            {
				int vertIndices[3];
                int edgeIndices[3];

                for (int j = 0; j < 3; j++)
                {
                    vertIndices[j] = mesh.faceVertex(i, (k + j) % 3);
                    edgeIndices[j] = mesh.faceEdge(i, (k + j) % 3);
                }

                if(mesh.edgeFace(edgeIndices[2], 0) == -1 || mesh.edgeFace(edgeIndices[2], 1) == -1)
					continue;

                for (int j = 0; j < edgeQuad.numQuadraturePoint(); j++)
                {
                    double s = edgeQuad.s(j);
                    double t = edgeQuad.t(j);
                    Eigen::Vector<double, 6> N = QuadraticElement::N(s, t);
                    Eigen::Matrix<double, 6, 2> dN = QuadraticElement::dN(s, t);

                    Eigen::Matrix<double, 18, 1> v;

                    for (int l = 0; l < 3; l++)
                    {
                        v.segment<3>(3 * l) = vertPos.row(vertIndices[l]).transpose();
                    }
                    for (int l = 0; l < 3; l++)
                    {
                        v.segment<3>(9 + 3 * l) = edgePos.row(edgeIndices[l]).transpose();
                    }

                    Eigen::Matrix<double, 6, 18> dxdr;
                    dxdr.setZero();
                    for (int m = 0; m < 3; m++)
                    {
                        for (int l = 0; l < 6; l++)
                        {
                            dxdr(m, 3 * l + m) += dN(l, 0);
                            dxdr(3 + m, 3 * l + m) += dN(l, 1);
                        }
                    }

                    Eigen::Vector3d pos(0, 0, 0);
                    for(int l=0; l<6; l++)
						pos += N(l) * v.segment<3>(3 * l);
					framePts.row(idx) = pos.transpose();

                    Eigen::Matrix<double, 6, 1> dr = dxdr * v;
					frameUs.row(idx) = dr.block<3, 1>(0, 0).transpose().normalized();
					frameVs.row(idx) = dr.block<3, 1>(3, 0).transpose().normalized();
                    Eigen::Vector3d normal = faceAreaNormal<18>(dr, dxdr, NULL, NULL);
					frameNs.row(idx) = normal.transpose().normalized();
                    idx++;
                }
			}
        }
    }

    double elasticEnergy(const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& paramDomain,
        const Eigen::MatrixXd* restPos,
        const Eigen::MatrixXd* restEdgePos,
        const Eigen::MatrixXd& curPos,
        const Eigen::MatrixXd& curEdgePos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        std::vector<bool> *violatedConstraints,
        Eigen::VectorXd* deriv,
        std::vector<Eigen::Triplet<double> >* hess
    )
    {
        double faceEnergy = 0;
		double edgeEnergy = 0;
        double barrierEnergy = 0;
        int nfaces = mesh.nFaces();
        int nverts = curPos.rows();
        int nedges = mesh.nEdges();

        if(violatedConstraints)
            violatedConstraints->clear();        

        if (deriv)
        {
            deriv->resize(3 * (nverts + nedges));
            deriv->setZero();
        }
        if (hess)
            hess->clear();
        
        
        {
            // Face energy 
            
            FourthOrderEdgeBasedTriangleQuadrature quad;
            
            for (int i = 0; i < nfaces; i++)
            {
                for (int j = 0; j < quad.numQuadraturePoint(); j++)
                {
                    double s = quad.s(j);
                    double t = quad.t(j);
					Eigen::Vector<double, 6> N = QuadraticElement::N(s, t);
					Eigen::Matrix<double, 6, 2> dN = QuadraticElement::dN(s, t);
					Eigen::Matrix<double, 6, 3> d2N = QuadraticElement::d2N(s, t);
                   
                    int vertIndices[3];
                    int edgeIndices[3];
                    Eigen::Matrix<double, 18, 1> v;
                    Eigen::Matrix<double, 18, 1> vbar;
                    for (int k = 0; k < 3; k++)
                    {
                        vertIndices[k] = mesh.faceVertex(i, k);
                        v.segment<3>(3 * k) = curPos.row(vertIndices[k]).transpose();
                        edgeIndices[k] = mesh.faceEdge(i, k);
                        v.segment<3>(9 + 3 * k) = curEdgePos.row(edgeIndices[k]).transpose();
                    }
                    if (restPos)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            vbar.segment<3>(3 * k) = restPos->row(mesh.faceVertex(i, k)).transpose();
                            vbar.segment<3>(9 + 3 * k) = restEdgePos->row(mesh.faceEdge(i, k)).transpose();
                        }
                    }
                    else
                    {
                        vbar = v;
                    }
                    

                    Eigen::Matrix<double, 6, 18> dxdr;
                    dxdr.setZero();
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 6; l++)
                        {
                            dxdr(k, 3 * l + k) = dN(l, 0);
                            dxdr(3 + k, 3 * l + k) = dN(l, 1);
                        }
                    }
                    std::vector<Eigen::Matrix<double, 3, 18> > dxd2r(4);
                    for (int k = 0; k < 4; k++)
                        dxd2r[k].setZero();
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 6; l++)
                        {
                            dxd2r[0](k, 3 * l + k) = d2N(l, 0);
                            dxd2r[1](k, 3 * l + k) = d2N(l, 1);
                            dxd2r[2](k, 3 * l + k) = d2N(l, 1);
                            dxd2r[3](k, 3 * l + k) = d2N(l, 2);
                        }
                    }

                    Eigen::Matrix<double, 6, 1> dr = dxdr * v;
                    Eigen::Matrix<double, 6, 1> drbar = dxdr * vbar;
                    std::vector<Eigen::Vector3d> d2r(4);
                    std::vector<Eigen::Vector3d> d2rbar(4);
                    for (int k = 0; k < 4; k++)
                    {
                        d2r[k] = dxd2r[k] * v;
                        d2rbar[k] = dxd2r[k] * vbar;
                    }

                    Eigen::Vector4d bbar;
                    if (restPos)
                    {
                        bbar = secondFundamentalForm<18>(drbar, d2rbar, dxdr, dxd2r, NULL, NULL);
                    }
                    else
                    {
                        bbar.setZero();
                    }                   

                    Eigen::Vector4d abar;
					Eigen::Vector3d parame0 = paramDomain.row(vertIndices[1]) - paramDomain.row(vertIndices[0]);
					Eigen::Vector3d parame1 = paramDomain.row(vertIndices[2]) - paramDomain.row(vertIndices[0]);

					abar << parame0.dot(parame0), parame0.dot(parame1), parame0.dot(parame1), parame1.dot(parame1);

                    double dA = 0.5 * std::sqrt(abar[0] * abar[3] - abar[1] * abar[2]);
                    
                    double coeff = quad.weight(j) * thickness * thickness * thickness / 12.0;

                    Eigen::Matrix<double, 1, 18> elementDeriv;
                    Eigen::Matrix<double, 18, 18> elementHess;

                    double h2norm =
                        MeanCurvatureSquaredEnergy<18>(dr, d2r, dxdr, dxd2r, abar, bbar, lameAlpha, lameBeta, deriv ? &elementDeriv : NULL, hess ? &elementHess : NULL);

                    faceEnergy += h2norm * coeff * dA;

                    if (deriv)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            deriv->segment<3>(3 * vertIndices[k]) += elementDeriv.segment<3>(3 * k) * coeff * dA;
                            deriv->segment<3>(3 * nverts + 3 * edgeIndices[k]) += elementDeriv.segment<3>(9 + 3 * k) * coeff * dA;   
                        }
                    }
                    if (hess)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                for (int m = 0; m < 3; m++)
                                {
                                    for (int n = 0; n < 3; n++)
                                    {
                                        hess->emplace_back(3 * vertIndices[k] + m, 3 * vertIndices[l] + n, elementHess(3 * k + m, 3 * l + n) * coeff * dA);
                                        hess->emplace_back(3 * vertIndices[k] + m, 3 * nverts + 3 * edgeIndices[l] + n, elementHess(3 * k + m, 9 + 3 * l + n) * coeff * dA);
                                        hess->emplace_back(3 * nverts + 3 * edgeIndices[k] + m, 3 * vertIndices[l] + n, elementHess(9 + 3 * k + m, 3 * l + n) * coeff * dA);
                                        hess->emplace_back(3 * nverts + 3 * edgeIndices[k] + m, 3 * nverts + 3 * edgeIndices[l] + n, elementHess(9 + 3 * k + m, 9 + 3 * l + n) * coeff * dA);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        
        {
            // Edge energy 
            
            SecondOrderEdgeQuadrature quad;
            
            for (int i = 0; i < nedges; i++)
            {
                if (mesh.edgeFace(i, 0) == -1 || mesh.edgeFace(i, 1) == -1)
                    continue;

                for (int j = 0; j < quad.numQuadraturePoint(); j++)
                {
                    double s1 = quad.s(j);
                    double s2 = quad.s(quad.numQuadraturePoint() - 1 - j);
                    double t1 = quad.t(j);
                    double t2 = quad.t(quad.numQuadraturePoint() - 1 - j);
					Eigen::Vector<double, 6> N = QuadraticElement::N(s1, t1);
					Eigen::Vector<double, 6> N2 = QuadraticElement::N(s2, t2);
					Eigen::Matrix<double, 6, 2> dN1 = QuadraticElement::dN(s1, t1);
					Eigen::Matrix<double, 6, 2> dN2 = QuadraticElement::dN(s2, t2);
                    
                    int vertIndices[4];
                    int edgeIndices[5];
                    Eigen::Matrix<double, 27, 1> v;
                    Eigen::Matrix<double, 27, 1> vbar;
                    int fv0 = mesh.edgeVertex(i, 0);
                    int fv1 = mesh.edgeVertex(i, 1);
                    int f0 = mesh.edgeFace(i, 0);
                    int f1 = mesh.edgeFace(i, 1);
                    vertIndices[0] = fv0;
                    vertIndices[1] = fv1;
                    vertIndices[2] = mesh.edgeOppositeVertex(i, 0);
                    vertIndices[3] = mesh.edgeOppositeVertex(i, 1);
                    edgeIndices[0] = i;
                    for (int k = 0; k < 3; k++)
                    {
                        if (mesh.faceVertex(f0, k) == fv0)
                            edgeIndices[1] = mesh.faceEdge(f0, k);
                        if(mesh.faceVertex(f0, k) == fv1)
							edgeIndices[2] = mesh.faceEdge(f0, k);
                        if(mesh.faceVertex(f1, k) == fv1)
                            edgeIndices[3] = mesh.faceEdge(f1, k);
						if (mesh.faceVertex(f1, k) == fv0)
							edgeIndices[4] = mesh.faceEdge(f1, k);
                    }

                    int localIndices1[6] = { 0, 1, 2, 5, 6, 4 };
                    int localIndices2[6] = { 1, 0, 3, 7, 8, 4 };

                    for (int k = 0; k < 4; k++)
                    {
                        v.segment<3>(3 * k) = curPos.row(vertIndices[k]).transpose();                        
                    }
                    for (int k = 0; k < 5; k++)
                    {
                        v.segment<3>(12 + 3 * k) = curEdgePos.row(edgeIndices[k]).transpose();
                    }
                    if (restPos)
                    {
                        for (int k = 0; k < 4; k++)
                        {
                            vbar.segment<3>(3 * k) = restPos->row(vertIndices[k]).transpose();                            
                        }
                        for (int k = 0; k < 5; k++)
                        {
                            vbar.segment<3>(12 + 3 * k) = restEdgePos->row(edgeIndices[k]).transpose();
                        }
                    }
                    else
                    {
                        vbar = v;
                    }
                    
                    Eigen::Matrix<double, 6, 27> dxdr1;
                    Eigen::Matrix<double, 6, 27> dxdr2;
                    dxdr1.setZero();
                    dxdr2.setZero();
                    for (int k = 0; k < 3; k++)
                    {
                        for (int l = 0; l < 6; l++)
                        {
                            dxdr1(k, 3 * localIndices1[l] + k) += dN1(l, 0);
                            dxdr1(3 + k, 3 * localIndices1[l] + k) += dN1(l, 1);
                            dxdr2(k, 3 * localIndices2[l] + k) += dN2(l, 0);
                            dxdr2(3 + k, 3 * localIndices2[l] + k) += dN2(l, 1);
                        }
                    }

                    double thetabar = 0;
                    if (restPos)
                    {
                        Eigen::Matrix<double, 6, 1> drbar1 = dxdr1 * vbar;
                        Eigen::Matrix<double, 6, 1> drbar2 = dxdr2 * vbar;

                        Eigen::Vector3d restNormal1 = faceAreaNormal<27>(drbar1, dxdr1, NULL, NULL);
                        Eigen::Vector3d restNormal2 = faceAreaNormal<27>(drbar2, dxdr2, NULL, NULL);

                        thetabar = LibShell::angle(restNormal1, restNormal2, drbar1.segment<3>(0), NULL, NULL);
                    }

                    Eigen::Matrix<double, 6, 1> dr1 = dxdr1 * v;
                    Eigen::Matrix<double, 6, 1> dr2 = dxdr2 * v;

                    Eigen::Matrix<double, 3, 27> nderiv1;
                    std::vector<Eigen::Matrix<double, 27, 27> > nhess1(3);
                    Eigen::Vector3d normal1 = faceAreaNormal<27>(dr1, dxdr1, &nderiv1, &nhess1);

                    Eigen::Matrix<double, 3, 27> nderiv2;
                    std::vector<Eigen::Matrix<double, 27, 27> > nhess2(3);
                    Eigen::Vector3d normal2 = faceAreaNormal<27>(dr2, dxdr2, &nderiv2, &nhess2);

                    Eigen::Vector3d axis = dr1.segment<3>(0);
                    Eigen::Matrix<double, 3, 27> axisderiv;
                    for(int k=0; k<3; k++)
                    {
                        axisderiv.row(k) = dxdr1.row(k);
                    }
                    
                    Eigen::Matrix<double, 1, 9> thetaderiv;
                    Eigen::Matrix<double, 9, 9> thetahess;
                    double theta = LibShell::angle(normal1, normal2, axis, &thetaderiv, &thetahess);

                    double weight = thickness * thickness * thickness / 12.0 * (lameAlpha + 2.0 * lameBeta);

					Eigen::Vector3d e0bar = paramDomain.row(fv1) - paramDomain.row(fv0);
					Eigen::Vector3d e1bar = paramDomain.row(vertIndices[2]) - paramDomain.row(fv0);
					Eigen::Vector3d e2bar = paramDomain.row(vertIndices[3]) - paramDomain.row(fv0);

                    double A1bar = e0bar.cross(e1bar).norm() / 2.0;
                    double A2bar = e0bar.cross(e2bar).norm() / 2.0;
                    
                    double coeff = 3.0 * quad.weight(j) * e0bar.squaredNorm() / (A1bar + A2bar) * weight;

                    double dscale1 = 1e4 / A1bar / A1bar;
                    double degen1 = 0.5 * normal1.squaredNorm() * dscale1;
                    double dscale2 = 1e4 / A2bar / A2bar;
                    double degen2 = 0.5 * normal2.squaredNorm() * dscale2;
                    double d1deriv = 0;
                    double d1hess = 0;
                    double d2deriv = 0;
                    double d2hess = 0;
                    double barrier1 = barrierPotential(degen1, d1deriv, d1hess);
                    if (violatedConstraints)
                        violatedConstraints->push_back(barrier1 > 0);
                    
					double barrier2 = barrierPotential(degen2, d2deriv, d2hess);
                    if (violatedConstraints)
                        violatedConstraints->push_back(barrier2 > 0);

                    barrierEnergy += barrier1 + barrier2;

                    edgeEnergy += 0.5 * coeff * (theta - thetabar) * (theta - thetabar);

                    if (deriv)
                    {
                        Eigen::Matrix<double, 1, 27> elementDeriv = thetaderiv.segment<3>(0) * nderiv1 + thetaderiv.segment<3>(3) * nderiv2 + thetaderiv.segment<3>(6) * axisderiv;
                        Eigen::Matrix<double, 1, 27> barrierDeriv = dscale1 * d1deriv * normal1.transpose() * nderiv1 + dscale2 * d2deriv * normal2.transpose() * nderiv2;
                        for (int k = 0; k < 4; k++)
                        {
                            deriv->segment<3>(3 * vertIndices[k]) += elementDeriv.segment<3>(3 * k) * coeff * (theta - thetabar);
                            deriv->segment<3>(3 * vertIndices[k]) += barrierDeriv.segment<3>(3 * k);
                        }
                        for(int k=0; k<5; k++)
                        {
                            deriv->segment<3>(3 * nverts + 3 * edgeIndices[k]) += elementDeriv.segment<3>(12 + 3 * k) * coeff * (theta - thetabar);                            
                            deriv->segment<3>(3 * nverts + 3 * edgeIndices[k]) += barrierDeriv.segment<3>(12 + 3 * k);
                        }
                    }
                    if (hess)
                    {
                        Eigen::Matrix<double, 27, 27> barrierHess = dscale1 * dscale1 * nderiv1.transpose() * normal1 * d1hess * normal1.transpose() * nderiv1;
                        barrierHess += dscale2 * dscale2 * nderiv2.transpose() * normal2 * d2hess * normal2.transpose() * nderiv2;
                        barrierHess += dscale1 * d1deriv * nderiv1.transpose() * nderiv1;
                        barrierHess += dscale2 * d2deriv * nderiv2.transpose() * nderiv2;
                        for (int k = 0; k < 3; k++)
                        {
                            barrierHess += dscale1 * d1deriv * normal1[k] * nhess1[k];
                            barrierHess += dscale2 * d2deriv * normal2[k] * nhess2[k];
                        }

                        Eigen::Matrix<double, 1, 27> elementDeriv = thetaderiv.segment<3>(0) * nderiv1 + thetaderiv.segment<3>(3) * nderiv2 + thetaderiv.segment<3>(6) * axisderiv;
						Eigen::Matrix<double, 27, 27> elementHess = elementDeriv.transpose() * elementDeriv;
						elementHess += (theta - thetabar) * nderiv1.transpose() * thetahess.block<3, 3>(0, 0) * nderiv1;
						elementHess += (theta - thetabar) * nderiv1.transpose() * thetahess.block<3, 3>(0, 3) * nderiv2;
                        elementHess += (theta - thetabar) * nderiv1.transpose() * thetahess.block<3, 3>(0, 6) * axisderiv;
                        elementHess += (theta - thetabar) * nderiv2.transpose() * thetahess.block<3, 3>(3, 0) * nderiv1;
                        elementHess += (theta - thetabar) * nderiv2.transpose() * thetahess.block<3, 3>(3, 3) * nderiv2;
                        elementHess += (theta - thetabar) * nderiv2.transpose() * thetahess.block<3, 3>(3, 6) * axisderiv;
                        elementHess += (theta - thetabar) * axisderiv.transpose() * thetahess.block<3, 3>(6, 0) * nderiv1;
                        elementHess += (theta - thetabar) * axisderiv.transpose() * thetahess.block<3, 3>(6, 3) * nderiv2;
                        elementHess += (theta - thetabar) * axisderiv.transpose() * thetahess.block<3, 3>(6, 6) * axisderiv;
												
                        for (int k = 0; k < 3; k++)
                        {
                            elementHess += (theta - thetabar) * thetaderiv[k] * nhess1[k];
                            elementHess += (theta - thetabar) * thetaderiv[3 + k] * nhess2[k];
                        }

                        for (int m = 0; m < 3; m++)
                        {
                            for (int n = 0; n < 3; n++)
                            {
                                for (int k = 0; k < 4; k++)
                                {
                                    for (int l = 0; l < 4; l++)
                                    {
                                        if (elementHess(3 * k + m, 3 * l + n) != 0)
                                        {
                                            hess->emplace_back(3 * vertIndices[k] + m, 3 * vertIndices[l] + n, elementHess(3 * k + m, 3 * l + n) * coeff);
                                        }
                                        if (barrierHess(3 * k + m, 3 * l + n) != 0)
                                        {
                                            hess->emplace_back(3 * vertIndices[k] + m, 3 * vertIndices[l] + n, barrierHess(3 * k + m, 3 * l + n));
                                        }
                                    }
                                    for (int l = 0; l < 5; l++)
                                    {
                                        if (elementHess(3 * k + m, 12 + 3 * l + n) != 0)
                                        {
                                            hess->emplace_back(3 * vertIndices[k] + m, 3 * nverts + 3 * edgeIndices[l] + n, elementHess(3 * k + m, 12 + 3 * l + n) * coeff);
                                        }
                                        if (barrierHess(3 * k + m, 12 + 3 * l + n) != 0)
                                        {
                                            hess->emplace_back(3 * vertIndices[k] + m, 3 * nverts + 3 * edgeIndices[l] + n, barrierHess(3 * k + m, 12 + 3 * l + n));
                                        }
                                    }
                                }
                                for (int k = 0; k < 5; k++)
                                {
                                    for (int l = 0; l < 4; l++)
                                    {
                                        if (elementHess(12 + 3 * k + m, 3 * l + n) != 0)
                                        {
                                            hess->emplace_back(3 * nverts + 3 * edgeIndices[k] + m, 3 * vertIndices[l] + n, elementHess(12 + 3 * k + m, 3 * l + n) * coeff);
                                        }
                                        if (barrierHess(12 + 3 * k + m, 3 * l + n) != 0)
                                        {
                                            hess->emplace_back(3 * nverts + 3 * edgeIndices[k] + m, 3 * vertIndices[l] + n, barrierHess(12 + 3 * k + m, 3 * l + n));
                                        }
                                    }
                                    for (int l = 0; l < 5; l++)
                                    {
                                        if (elementHess(12 + 3 * k + m, 12 + 3 * l + n) != 0)
                                        {
                                            hess->emplace_back(3 * nverts + 3 * edgeIndices[k] + m, 3 * nverts + 3 * edgeIndices[l] + n, elementHess(12 + 3 * k + m, 12 + 3 * l + n) * coeff);
                                        }
                                        if (barrierHess(12 + 3 * k + m, 12 + 3 * l + n) != 0)
                                        {
                                            hess->emplace_back(3 * nverts + 3 * edgeIndices[k] + m, 3 * nverts + 3 * edgeIndices[l] + n, barrierHess(12 + 3 * k + m, 12 + 3 * l + n));
                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }
		std::cout << "Face energy: " << faceEnergy << ", edge energy: " << edgeEnergy << ", barrier energy: " << barrierEnergy << std::endl;
        return faceEnergy + edgeEnergy + barrierEnergy;
    }

    void optimizeMidedgeDOFs(const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& paramDomain,
        const Eigen::MatrixXd* restPos,
        const Eigen::MatrixXd* restEdgePos,
        const Eigen::MatrixXd& curPos,
        Eigen::MatrixXd& edgePosGuess,
        double thickness,
        double lameAlpha,
        double lameBeta)
    {
		int nverts = curPos.rows();
		int nedges = mesh.nEdges();
        
        std::vector<Eigen::Triplet<double> > Pcoeffs;
        int reducedDofs = 0;
        for (int i = 0; i < nedges; i++)
        {
            //if(mesh.edgeFace(i,0) == -1 || mesh.edgeFace(i,1) == -1)
				//continue;
            for (int j = 0; j < 3; j++)
            {
                Pcoeffs.push_back(Eigen::Triplet<double>(reducedDofs, 3 * nverts + 3 * i + j, 1.0));
                reducedDofs++;
            }
        }
		Eigen::SparseMatrix<double> P(reducedDofs, 3 * (nverts + nedges));
		P.setFromTriplets(Pcoeffs.begin(), Pcoeffs.end());

        double tol = 1e-5;
        
        std::vector<Eigen::Triplet<double> > Icoeffs;
        for (int i = 0; i < reducedDofs; i++)
        {
            Icoeffs.push_back({ i, i, 1.0 });
        }
        Eigen::SparseMatrix<double> I(reducedDofs, reducedDofs);
        I.setFromTriplets(Icoeffs.begin(), Icoeffs.end());

        Eigen::SparseMatrix<double> PT = P.transpose();

        double reg = 1e-6;
        while (true)
        {
            std::vector<Eigen::Triplet<double> > Hcoeffs;
            Eigen::VectorXd F;
            double origEnergy = elasticEnergy(mesh, paramDomain, restPos, restEdgePos, curPos, edgePosGuess, thickness, lameAlpha, lameBeta, NULL, &F, &Hcoeffs);
            Eigen::VectorXd PF = P * F;
            std::cout << "Force resid now: " << PF.norm() << ", energy: " << origEnergy << ", reg: " << reg << std::endl;
            if (PF.norm() < tol)
                return;
            Eigen::SparseMatrix<double> H(3 * nverts + 3 * nedges, 3 * nverts + 3 * nedges);
            H.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());
            Eigen::SparseMatrix<double> PHPT = P * H * PT;
            Eigen::SparseMatrix<double> M = PHPT + reg * I;
            Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver(M);
            Eigen::VectorXd update = solver.solve(-PF);
            if (solver.info() != Eigen::Success) {
                std::cout << "Solve failed" << std::endl;
                reg *= 2.0;
                continue;
            }
            Eigen::VectorXd fullUpdate = PT * update;
            Eigen::MatrixXd newedgeDOFs = edgePosGuess;
            for (int i = 0; i < nedges; i++)
            {
                newedgeDOFs.row(i) += fullUpdate.segment<3>(3 * nverts + 3 * i).transpose();
            }
            double newenergy = elasticEnergy(mesh, paramDomain, restPos, restEdgePos, curPos, newedgeDOFs, thickness, lameAlpha, lameBeta, NULL, NULL, NULL);
            if (newenergy > origEnergy)
            {
                std::cout << "Not a descent step, " << origEnergy << " -> " << newenergy << std::endl;
                reg *= 2.0;
                continue;
            }

            edgePosGuess = newedgeDOFs;
            reg *= 0.5;
        }
    }

    double optimalElasticEnergy(
        const LibShell::MeshConnectivity& mesh,        
        const Eigen::MatrixXd* restPos,
        const Eigen::MatrixXd* restEdgePos,
        const Eigen::MatrixXd& curPos,
        const Eigen::MatrixXd& curEdgePos,
        double thickness,
        double lameAlpha,
        double lameBeta)
    {
        Eigen::MatrixXd paramDomain;
        if (restPos)
            paramDomain = *restPos;
        else
            paramDomain = curPos;

		Eigen::MatrixXd restEdgePosGuess;
        if (restPos)
        {
            restEdgePosGuess = *restEdgePos;
			optimizeMidedgeDOFs(mesh, paramDomain, NULL, NULL, *restPos, restEdgePosGuess, thickness, lameAlpha, lameBeta);
        }

		Eigen::MatrixXd curEdgePosGuess = curEdgePos;        
        optimizeMidedgeDOFs(mesh, paramDomain, restPos, restPos ? &restEdgePosGuess : NULL, curPos, curEdgePosGuess, thickness, lameAlpha, lameBeta);               

        /*
        Eigen::MatrixXd visV;
        Eigen::MatrixXi visF;
        Eigen::VectorXd visEnergyDensity;
        Eigen::MatrixXd framePts;
        Eigen::MatrixXd frameUs;
        Eigen::MatrixXd frameVs;
        Eigen::MatrixXd frameNs;
        visualize(mesh, paramDomain, curPos, curEdgePosGuess, thickness, lameAlpha, lameBeta, visV, visF, visEnergyDensity, framePts, frameUs, frameVs, frameNs);
        writeQuadElements("failed.dat", curPos, curEdgePosGuess, visV, visF, visEnergyDensity, framePts, frameNs);
        */

        double result = elasticEnergy(mesh, paramDomain, restPos, restPos ? &restEdgePosGuess : NULL, curPos, curEdgePosGuess, thickness, lameAlpha, lameBeta, NULL, NULL, NULL);        
        return result;
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
        return optimalElasticEnergy(mesh, &restPos, &restEdgePos, curPos, curEdgePos, thickness, lameAlpha, lameBeta);
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
        return optimalElasticEnergy(mesh, &restPos, &restEdgePos, curPos, curEdgePos, thickness, lameAlpha, lameBeta);
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
        return optimalElasticEnergy(mesh, &restPos, &restEdgePos, curPos, curEdgePos, thickness, lameAlpha, lameBeta);
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
        return optimalElasticEnergy(mesh, NULL, NULL, curPos, curEdgePos, thickness, lameAlpha, lameBeta);
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
        return optimalElasticEnergy(mesh, NULL, NULL, curPos, curEdgePos, thickness, lameAlpha, lameBeta);
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
        return optimalElasticEnergy(mesh, &restPos, &restEdgePos, curPos, curEdgePos, thickness, lameAlpha, lameBeta);
    }
};

#endif