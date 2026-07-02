#ifndef DKTBENDING_H
#define DKTBENDING_H

#include "Model.h"
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>

static Eigen::Vector3d perpTo(const Eigen::Vector3d& v)
{
    int smallestidx = 0;
    for (int i = 1; i < 3; i++)
    {
        if(std::fabs(v[i]) < std::fabs(v[smallestidx]))
        {
            smallestidx = i;
		}
    }
    Eigen::Vector3d axis(0, 0, 0);
    axis[smallestidx] = 1.0;
	return v.cross(axis).normalized();
}

class DKTBending : public Model
{
private:
    Eigen::Matrix<double, 3, 9> secondFundamentalForm(const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& curPos,
        int faceIdx,
        double s,
        double t
    )
    {
        Eigen::Vector<double, 6> N;
        N << 2.0 * (1.0 - s - t) * (0.5 - s - t),
            s* (2.0 * s - 1.0),
            t* (2.0 * t - 1.0),
            4.0 * s * t,
            4.0 * t * (1.0 - s - t),
            4.0 * s * (1.0 - s - t);
        Eigen::Matrix<double, 6, 2> dN;
        dN << -3.0 + 4.0 * s + 4.0 * t, -3.0 + 4.0 * s + 4.0 * t,
            4.0 * s - 1.0, 0.0,
            0, 4.0 * t - 1.0,
            4.0 * t, 4.0 * s,
            -4.0 * t, 4.0 * (1.0 - s - 2.0 * t),
            4.0 * (1.0 - 2.0 * s - t), -4.0 * s;

        Eigen::Vector3d v[3];
        for (int i = 0; i < 3; i++)
        {
            v[i] = curPos.row(mesh.faceVertex(faceIdx, i)).transpose();
        }
        Eigen::Vector3d e[3];
        for (int i = 0; i < 3; i++)
        {
            e[i] = v[(i + 2) % 3] - v[(i + 1) % 3];
        }
        Eigen::Vector3d faceN = e[0].cross(e[1]);
        faceN.normalize();

        Eigen::Vector3d tan[3];
        Eigen::Vector3d n[3];
        for (int i = 0; i < 3; i++)
        {
            tan[i] = e[i].normalized();
            n[i] = faceN.cross(tan[i]);
        }

        Eigen::Matrix<double, 18, 9> P;
        P.setZero();
        for (int i = 0; i < 3; i++)
        {
            P.block<3, 3>(3 * i, 3 * i) = Eigen::Matrix3d::Identity();
        }
        for (int i = 0; i < 3; i++)
        {
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            P.block<3, 3>(9 + 3 * i, 3 * j) = 0.5 * n[i] * n[i].transpose() - 0.25 * tan[i] * tan[i].transpose();
            P.block<3, 3>(9 + 3 * i, 3 * k) = 0.5 * n[i] * n[i].transpose() - 0.25 * tan[i] * tan[i].transpose();
        }
        Eigen::Matrix<double, 3, 9> bentries;
        bentries.setZero();
        Eigen::Vector3d dr1 = v[1] - v[0];
        Eigen::Vector3d dr2 = v[2] - v[0];
        for (int i = 0; i < 6; i++)
        {
            bentries.row(0) += dN(i, 0) * dr1.transpose() * P.block<3, 9>(3 * i, 0);
            bentries.row(1) += 0.5 * dN(i, 1) * dr1.transpose() * P.block<3, 9>(3 * i, 0);
            bentries.row(1) += 0.5 * dN(i, 0) * dr2.transpose() * P.block<3, 9>(3 * i, 0);
            bentries.row(2) += dN(i, 1) * dr2.transpose() * P.block<3, 9>(3 * i, 0);
        }
        return bentries;
    }

    void localBendingEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd* restPos,
        const Eigen::VectorXd* restDirectors,
        const Eigen::MatrixXd& curPos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        int faceIdx,
        Eigen::Matrix<double, 9, 9>& quadraticPart,
        Eigen::Vector<double, 9>& linearPart,
        double& constantPart)
    {
        quadraticPart.setZero();
        linearPart.setZero();
        constantPart = 0;

        Eigen::Vector3d abarentriesVec;
        if (restPos)
        {
            Eigen::Vector3d v[3];
            for (int i = 0; i < 3; i++)
            {
                v[i] = restPos->row(mesh.faceVertex(faceIdx, i)).transpose();
            }
            abarentriesVec[0] = (v[1] - v[0]).dot(v[1] - v[0]);
            abarentriesVec[1] = (v[1] - v[0]).dot(v[2] - v[0]);
            abarentriesVec[2] = (v[2] - v[0]).dot(v[2] - v[0]);
        }
        else
        {
            Eigen::Vector3d v[3];
            for (int i = 0; i < 3; i++)
            {
                v[i] = curPos.row(mesh.faceVertex(faceIdx, i)).transpose();
            }
            abarentriesVec[0] = (v[1] - v[0]).dot(v[1] - v[0]);
            abarentriesVec[1] = (v[1] - v[0]).dot(v[2] - v[0]);
            abarentriesVec[2] = (v[2] - v[0]).dot(v[2] - v[0]);
        }

        // three-point quadrature
        double s[3] = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
        double t[3] = { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 };
        double weights[3] = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
        for (int i = 0; i < 3; i++)
        {
            Eigen::Matrix<double, 3, 9> bentries = secondFundamentalForm(mesh, curPos, faceIdx, s[i], t[i]);
            Eigen::Vector3d bbarentriesVec;
            if (restPos)
            {
                Eigen::Matrix<double, 3, 9> bbarentries = secondFundamentalForm(mesh, *restPos, faceIdx, s[i], t[i]);
                Eigen::Vector<double, 9> localDirectors;

                for (int j = 0; j < 3; j++)
                {
                    localDirectors.segment<3>(3 * j) = restDirectors->segment<3>(3 * mesh.faceVertex(faceIdx, j)).transpose();
                }
                bbarentriesVec = bbarentries * localDirectors;
            }
            else
            {
                bbarentriesVec.setZero();
            }

            double dA = 0.5 * std::sqrt(abarentriesVec[0] * abarentriesVec[2] - abarentriesVec[1] * abarentriesVec[1]);
            double coeff = weights[i] * thickness * thickness * thickness / 12.0 * dA;

            double a0 = abarentriesVec[0];
            double a1 = abarentriesVec[1];
            double a2 = abarentriesVec[2];
            double bb0 = bbarentriesVec[0];
            double bb1 = bbarentriesVec[1];
            double bb2 = bbarentriesVec[2];            
            
            double denom = (a1 * a1 - a0 * a2) * (a1 * a1 - a0 * a2);

            constantPart += coeff / 2.0 / denom *
                ((a2 * bb0 - 2.0 * a1 * bb1 + a0 * bb2) * (a2 * bb0 - 2.0 * a1 * bb1 + a0 * bb2) * lameAlpha 
                    + 2.0 * (a2 * a2 * bb0 * bb0 + 2.0 * a2 * bb1 * (a0 * bb1 - 2.0 * a1 * bb0) - 4.0 * a0 * a1 * bb1 * bb2 + a0 * a0 * bb2 * bb2 + 2.0 * a1 * a1 * (bb1 * bb1 + bb0 * bb2)) * lameBeta);

            Eigen::Vector3d linearCoeffs;
            linearCoeffs <<
                -a0 * a2 * bb2 * lameAlpha - 2.0 * a1 * a1 * bb2 * lameBeta - a2 * a2 * bb0 * (lameAlpha + 2.0 * lameBeta) + 2.0 * a1 * a2 * bb1 * (lameAlpha + 2.0 * lameBeta),
                -4.0 * a0 * a2 * bb1 * lameBeta - 4.0 * a1 * a1 * bb1 * (lameAlpha + lameBeta) + 2.0 * a1 * (a2 * bb0 + a0 * bb2) * (lameAlpha + 2.0 * lameBeta),
                -a0 * a2 * bb0 * lameAlpha - 2.0 * a1 * a1 * bb0 * lameBeta - a0 * a0 * bb2 * (lameAlpha + 2.0 * lameBeta) + 2.0 * a0 * a1 * bb1 * (lameAlpha + 2.0 * lameBeta);
            linearPart += coeff / denom * linearCoeffs.transpose() * bentries;

            Eigen::Matrix3d quadraticCoeffs;
            quadraticCoeffs <<
                a2 * a2 * (lameAlpha + 2.0 * lameBeta), -2.0 * a1 * a2 * (lameAlpha + 2.0 * lameBeta), a0* a2* lameAlpha + 2.0 * a1 * a1 * lameBeta,
				-2.0 * a1 * a2 * (lameAlpha + 2.0 * lameBeta), 4.0 * (a0 * a2 * lameBeta + a1 * a1 * (lameAlpha + lameBeta)), -2.0 * a0 * a1 * (lameAlpha + 2.0 * lameBeta),                
                a0* a2* lameAlpha + 2.0 * a1 * a1 * lameBeta, -2.0 * a0 * a1 * (lameAlpha + 2.0 * lameBeta), a0* a0* (lameAlpha + 2.0 * lameBeta);            
            quadraticPart += coeff / 2.0 / denom * bentries.transpose() * quadraticCoeffs * bentries;
        }
    }

    void bendingEnergy(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd* restPos,
        const Eigen::VectorXd* restDirectors,
        const Eigen::MatrixXd& curPos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        Eigen::SparseMatrix<double>& quadraticPart,
        Eigen::VectorXd& linearPart,
        double& constantPart)
    {
        int nverts = curPos.rows();
        int nfaces = mesh.nFaces();
        quadraticPart.resize(3 * nverts, 3 * nverts);
        linearPart.resize(3 * nverts);
        linearPart.setZero();
        constantPart = 0;
        std::vector<Eigen::Triplet<double>> quadraticPartTriplets;
        for (int i = 0; i < nfaces; i++)
        {
            Eigen::Matrix<double, 9, 9> localQuadraticPart;
            Eigen::Vector<double, 9> localLinearPart;
            double localConstantPart;
            localBendingEnergy(mesh, restPos, restDirectors, curPos, thickness, lameAlpha, lameBeta,
                i, localQuadraticPart, localLinearPart, localConstantPart);
			double check = (localQuadraticPart - localQuadraticPart.transpose()).norm();
            if(check > 1e-6)
            {
                std::cout << "Local quadratic part is not symmetric, norm: " << check << std::endl;
                exit(-1);
			}
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            quadraticPartTriplets.push_back(
                                { 3 * mesh.faceVertex(i, j) + k, 3 * mesh.faceVertex(i, l) + m, localQuadraticPart(3 * j + k, 3 * l + m) });
                        }
                    }
                }
                linearPart.segment<3>(3 * mesh.faceVertex(i, j)) += localLinearPart.segment<3>(3 * j);
            }
            constantPart += localConstantPart;
        }
		quadraticPart.setFromTriplets(quadraticPartTriplets.begin(), quadraticPartTriplets.end());
    }

    void initializeDirectors(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd& curPos,
        Eigen::VectorXd& directors)
    {
		int nverts = curPos.rows();
        directors.resize(3 * nverts);
        directors.setZero();
		int nfaces = mesh.nFaces();
        for (int i = 0; i < nfaces; i++)
        {
            Eigen::Vector3d v[3];
            for (int j = 0; j < 3; j++)
            {
                v[j] = curPos.row(mesh.faceVertex(i, j)).transpose();
            }
            Eigen::Vector3d e1 = v[1] - v[0];
            Eigen::Vector3d e2 = v[2] - v[0];
            Eigen::Vector3d faceN = e1.cross(e2);
            faceN.normalize();
            for (int j = 0; j < 3; j++)
            {
                directors.segment<3>(3 * mesh.faceVertex(i, j)) += faceN;
			}
        }
        for (int i = 0; i < nverts; i++)
        {
            Eigen::Vector3d dir = directors.segment<3>(3 * i);
            if (dir.norm() < 1e-6)
            {
                dir.setRandom();
            }
            directors.segment<3>(3 * i) = dir.normalized();
		}
    }

    void optimizeDirectors(
        const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd* restPos,
        const Eigen::VectorXd* restDirectors,
        const Eigen::MatrixXd& curPos,
        double thickness,
        double lameAlpha,
        double lameBeta,
        Eigen::VectorXd &directorsGuess)
    {
		Eigen::SparseMatrix<double> quadraticPart;
        Eigen::VectorXd linearPart;
        double constantPart;
        bendingEnergy(mesh, restPos, restDirectors, curPos, thickness, lameAlpha, lameBeta,
			quadraticPart, linearPart, constantPart);

        int nverts = curPos.rows();
        assert(directorsGuess.size() == 3 * nverts);

        while (true)
        {
			std::cout << "Energy is: " << directorsGuess.transpose() * quadraticPart * directorsGuess + linearPart.dot(directorsGuess) + constantPart;
            std::cout << " = " << directorsGuess.transpose() * quadraticPart * directorsGuess << " + " 
				<< linearPart.dot(directorsGuess) << " + " << constantPart << std::endl;    
            for(int i=0; i<nverts; i++)
            {
                directorsGuess.segment<3>(3 * i).normalize();
			}
            std::vector<Eigen::Triplet<double>> Tcoeffs;
            for (int i = 0; i < nverts; i++)
            {
				Eigen::Vector3d dir = directorsGuess.segment<3>(3 * i);
				Eigen::Vector3d t1 = perpTo(dir);
				Eigen::Vector3d t2 = dir.cross(t1);
                for (int j = 0; j < 3; j++)
                {
					Tcoeffs.push_back({ 3 * i + j, 2 * i, t1[j] });
                    Tcoeffs.push_back({ 3 * i + j, 2 * i + 1, t2[j] });
                }
            }
			Eigen::SparseMatrix<double> T(3 * nverts, 2 * nverts);
            T.setFromTriplets(Tcoeffs.begin(), Tcoeffs.end());
            Eigen::SparseMatrix<double> H = T.transpose() * quadraticPart * T;
            Eigen::VectorXd b = -T.transpose() * quadraticPart * directorsGuess - 0.5 * T.transpose() * linearPart;
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(H);
            if (solver.info() != Eigen::Success)
            {
                std::cout << "Solve failed" << std::endl;
                exit(-1);
            }
            Eigen::VectorXd update = solver.solve(b);
            if (solver.info() != Eigen::Success)
            {
                std::cout << "Solve failed" << std::endl;
                exit(-1);
			}
            if (update.norm() < 1e-6)
            {
                break;
            }
			std::cout << "Update norm: " << update.norm() << std::endl;
            directorsGuess += T * update;
        }
		
    }

    double elasticEnergy(const LibShell::MeshConnectivity& mesh,
        const Eigen::MatrixXd *restPos,
        const Eigen::MatrixXd& curPos,
        double thickness,
        double lameAlpha,
        double lameBeta)
    {
        int nverts = curPos.rows();
        Eigen::VectorXd restDirectors(3 * nverts);
        if (restPos)
        {
            initializeDirectors(mesh, *restPos, restDirectors);
            optimizeDirectors(mesh, NULL, NULL, *restPos, thickness, lameAlpha, lameBeta, restDirectors);
            std::cout << "Done initializing rest directors" << std::endl;
        }
		Eigen::VectorXd directors(3 * nverts);
		initializeDirectors(mesh, curPos, directors);
		optimizeDirectors(mesh, restPos, restPos ? &restDirectors : NULL, curPos, thickness, lameAlpha, lameBeta, directors);
        std::cout << "Done initializing current directors" << std::endl;

		Eigen::SparseMatrix<double> quadraticPart;
        Eigen::VectorXd linearPart;
        double constantPart;
        bendingEnergy(mesh, restPos, restPos ? &restDirectors : NULL, curPos, thickness, lameAlpha, lameBeta,
            quadraticPart, linearPart, constantPart);
		return directors.transpose() * quadraticPart * directors + linearPart.dot(directors) + constantPart;
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