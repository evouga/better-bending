#include <Eigen/Core>
#include <Eigen/Sparse>

#include <igl/triangle/triangulate.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <sstream>
#include <iomanip>
#include <fstream>

#include "Model.h"
#include "BAC.h"
#include "DCS.h"
#include "QBCR.h"
#include "QBPL.h"
#include "DStheta.h"

void makeRectangularMesh(double W, double L, int Lverts, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
	int Wverts = std::max(2, int(W * Lverts / L));

    double triangleArea = 0.5 * W * L / (Lverts * Wverts);

    Eigen::MatrixXd Vin(2 * Wverts + 2 * Lverts, 2);
    Eigen::MatrixXi E(2 * Wverts + 2 * Lverts, 2);
    Eigen::MatrixXd dummyH(0, 2);
    Eigen::MatrixXd V2;
    Eigen::MatrixXi F2;

    int vrow = 0;
    int erow = 0;
    // top boundary
    for (int i = 1; i < Wverts; i++)
    {
        Vin(vrow, 0) = double(i) / double(Wverts) * W;
        Vin(vrow, 1) = L;
        if (i > 1)
        {
            E(erow, 0) = vrow - 1;
            E(erow, 1) = vrow;
            erow++;
        }
        vrow++;
    }
    // bottom boundary
    for (int i = 1; i < Wverts; i++)
    {
        Vin(vrow, 0) = double(i) / double(Wverts) * W;
        Vin(vrow, 1) = 0;
        if (i > 1)
        {
            E(erow, 0) = vrow - 1;
            E(erow, 1) = vrow;
            erow++;
        }
        vrow++;
    }
    // left boundary    
    for (int i = 0; i <= Lverts; i++)
    {
        Vin(vrow, 0) = 0;
        Vin(vrow, 1) = double(i) / double(Lverts) * L;
        if (i > 0)
        {
            E(erow, 0) = vrow - 1;
            E(erow, 1) = vrow;
            erow++;
        }
        vrow++;
    }
    // right boundary    
    for (int i = 0; i <= Lverts; i++)
    {
        Vin(vrow, 0) = W;
        Vin(vrow, 1) = double(i) / double(Lverts) * L;
        if (i > 0)
        {
            E(erow, 0) = vrow - 1;
            E(erow, 1) = vrow;
            erow++;
        }
        vrow++;
    }
    // missing four edges
    E(erow, 0) = (Wverts - 1) - 1;
    E(erow, 1) = 2 * (Wverts - 1) + 2 * (Lverts + 1) - 1;
    erow++;
    E(erow, 0) = 2 * (Wverts - 1) + (Lverts + 1);
    E(erow, 1) = 2 * (Wverts - 1) - 1;
    erow++;
    E(erow, 0) = Wverts - 1;
    E(erow, 1) = 2 * (Wverts - 1);
    erow++;
    E(erow, 0) = 2 * (Wverts - 1) + (Lverts + 1) - 1;
    E(erow, 1) = 0;
    erow++;

    assert(vrow == 2 * Lverts + 2 * Wverts);
    assert(erow == 2 * Lverts + 2 * Wverts);
    std::stringstream ss;
    ss << "a" << std::setprecision(30) << std::fixed << triangleArea << "qDY";
    igl::triangle::triangulate(Vin, E, dummyH, ss.str(), V2, F2);

    // roll up

    int nverts = V2.rows();

    V.resize(nverts, 3);
    
    for (int i = 0; i < nverts; i++)
    {
        Eigen::Vector2d q = V2.row(i).transpose();
        V(i, 0) = q[0];
        V(i, 1) = q[1];
        V(i, 2) = 0;        
    }
    
    F = F2;
}

void findBoundary(const Eigen::MatrixXd& V, double tol, std::vector<int>& pinned)
{
	int nverts = V.rows();
    for (int i = 0; i < nverts; i++)
    {
        if (V(i, 1) < tol)
            pinned.push_back(i);
    }
}

void computeVertexAreas(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXd& vertexAreas)
{
    int nverts = V.rows();
    vertexAreas.resize(nverts);
    vertexAreas.setZero();

    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d v0 = V.row(F(i, 0)).transpose();
        Eigen::Vector3d v1 = V.row(F(i, 1)).transpose();
        Eigen::Vector3d v2 = V.row(F(i, 2)).transpose();
        double area = 0.5 * (v1 - v0).cross(v2 - v0).norm();
        for (int j = 0; j < 3; j++)
        {
            vertexAreas[F(i, j)] += area / 3.0;
        }
    }
}

void optimizeDOFs(Model &model,
    Eigen::MatrixXd& curPos,
    Eigen::VectorXd& edgeDOFs,
    const std::vector<int> &pinned,
    const Eigen::VectorXd &extF,
    double tol)
{
    int nverts = curPos.rows();
    int nposdofs = nverts * 3;
    int nedgedofs = edgeDOFs.size();

    std::set<int> pinnedset;
    for (auto it : pinned)
        pinnedset.insert(it);

    std::vector<Eigen::Triplet<double> > Pcoeffs;

    int idx = 0;
    for (int i = 0; i < nverts; i++)
    {
        if (pinnedset.count(i))
            continue;
        for (int j = 0; j < 3; j++)
        {
            Pcoeffs.push_back({ idx, 3 * i + j, 1.0 });
            idx++;
        }
    }

    for (int i = 0; i < nedgedofs; i++)
    {
        Pcoeffs.push_back({ idx, nposdofs + i, 1.0 });
        idx++;
    }
    
    int reduceddofs = idx;

    std::vector<Eigen::Triplet<double> > Icoeffs;
    for (int i = 0; i < reduceddofs; i++)
    {        
        Icoeffs.push_back({ i, i, 1.0 });
    }
    Eigen::SparseMatrix<double> P(reduceddofs, nposdofs + nedgedofs);
    P.setFromTriplets(Pcoeffs.begin(), Pcoeffs.end());
    Eigen::SparseMatrix<double> I(reduceddofs, reduceddofs);
    I.setFromTriplets(Icoeffs.begin(), Icoeffs.end());

    Eigen::SparseMatrix<double> PT = P.transpose();


    double reg = 1e-6;
    while (true)
    {
        Eigen::VectorXd F;
        Eigen::SparseMatrix<double> H(nposdofs + nedgedofs, nposdofs + nedgedofs);

		double origEnergy = model.energy(curPos, edgeDOFs, &F, &H);
        origEnergy += -curPos.reshaped<Eigen::RowMajor>().dot(extF);

        F.segment(0, nposdofs) -= extF;
        Eigen::VectorXd PF = P * F;
        std::cout << "Force resid now: " << PF.norm() << ", energy: " << origEnergy << ", reg: " << reg;
        if (PF.norm() < tol)
        {
            std::cout << "; converged" << std::endl;
            return;
        }       
        
        Eigen::SparseMatrix<double> PHPT = P * H * PT;

        Eigen::SparseMatrix<double> M = PHPT + reg * I;
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver(M);
        Eigen::VectorXd update = solver.solve(-PF);
        if (solver.info() != Eigen::Success) {
            std::cout << "; solve failed" << std::endl;
            reg *= 2.0;
            continue;
        }
        else
        {
            std::cout << ", Newton decrement: " << update.norm() << std::endl;
        }

        if (update.norm() < 1e-12)
        {
            std::cout << "Stalled" << std::endl;
            return;
        }

        Eigen::VectorXd fullUpdate = PT * update;
        Eigen::MatrixXd newPos = curPos;
        for (int i = 0; i < nverts; i++)
        {
            newPos.row(i) += fullUpdate.segment<3>(3 * i).transpose();
        }

        Eigen::VectorXd newedgeDOFs = edgeDOFs + fullUpdate.segment(nposdofs, nedgedofs);

        if (!model.dofsValid(newPos, newedgeDOFs))
        {
            std::cout << "Edge DOFs invalid" << std::endl;
            reg *= 2.0;
            continue;
        }

		double newenergy = model.energy(newPos, newedgeDOFs, NULL, NULL);
        newenergy += -newPos.reshaped<Eigen::RowMajor>().dot(extF);

        if (newenergy > origEnergy)
        {
            std::cout << "Not a descent step, " << origEnergy << " -> " << newenergy << std::endl;
            reg *= 2.0;
            continue;
        }
        curPos = newPos;
        edgeDOFs = newedgeDOFs;
        reg *= 0.5;
    }
}

int main()
{
    std::string logfile = "log.txt";
    std::ofstream log(logfile);

//    polyscope::init();
    
    double nu = 0.35;
	double L = 1.0;
	double D = 1.0;
	double h = 1e-3;

    double density = 1; // doesn't matter, will divide out

	// D = young * thickness^3 / (12 * (1 - poisson^2))

	double Y = D * 12 * (1 - nu * nu) / (h * h * h);

    double lameAlpha = Y * nu / (1.0 - nu * nu);
    double lameBeta = Y / 2.0 / (1.0 + nu);

    double newtonTol = 1e-7;


	double initialW = 0.1;
	double dW = 0.1;
	double Wsteps = 9;

	for (int i = 0; i < Wsteps; i++)
	{
		double W = initialW + i * dW;

        log << "Testing W = " << W << std::endl;

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		makeRectangularMesh(W, L, 100, V, F);
        int nverts = V.rows();

        std::vector<int> pinned;
        findBoundary(V, 1.5 * L/100.0, pinned);

        Eigen::VectorXd vertexAreas;
        computeVertexAreas(V, F, vertexAreas);
        for (auto it : pinned)
            vertexAreas[it] = 0;

        std::map<std::string, Model*> models;
        models["BAC"] = new BAC(V, F, h, lameAlpha, lameBeta);
		models["DCS"] = new DCS(V, F, h, lameAlpha, lameBeta);
        models["QBPL"] = new QBPL(V, F, h, lameAlpha, lameBeta);
		models["QBCR"] = new QBCR(V, F, h, lameAlpha, lameBeta);
		models["DS(theta)"] = new DSTheta(V, F, h, lameAlpha, lameBeta);
        
        Eigen::MatrixXd curV = V;
        
        Eigen::VectorXd perturbation(nverts);
        perturbation.setRandom();
        for (auto it : pinned)
            perturbation[it] = 0;
        for (int i = 0; i < nverts; i++)
        {
            curV(i, 2) += 1e-4 * std::fabs(perturbation[i]);
        }

		std::map<std::string, Eigen::MatrixXd> curPoss;
		std::map<std::string, Eigen::VectorXd> curEdgeDOFs;

        log << "Models: ";
        bool first = true;
		for (auto it : models)
        {
            curPoss[it.first] = curV;
            curEdgeDOFs[it.first] = it.second->restEdgeDOFs();

            if (!first)
                log << ", ";
            else
                first = false;

			log << it.first;
        }
        log << std::endl;
        
        for (double gammaStar = 40; gammaStar >= 10; gammaStar -= 0.5)
        {
            double g = D * gammaStar / (density * L * L * L * h);

            Eigen::VectorXd extF(3 * nverts);
            extF.setZero();
            for (int i = 0; i < nverts; i++)
            {
                extF(3 * i + 0) = g * h * vertexAreas[i] * density;
            }

            std::cout << "==========" << std::endl;
            std::cout << "Starting solve: w/L = " << W << ", gammaStar = " << gammaStar << std::endl;
            std::cout << "==========" << std::endl;
            log << gammaStar << ": ";
            bool first = true;
            for (auto it : models)
            {
                if (!first)
                    log << ", ";
                else
                    first = false;

				optimizeDOFs(*it.second, curPoss[it.first], curEdgeDOFs[it.first], pinned, extF, newtonTol);
                double maxLateral = 0;
                for (int i = 0; i < nverts; i++)
                {
                    maxLateral = std::max(maxLateral, std::fabs(curPoss[it.first](i, 2)));
				}
                log << maxLateral;
            }
            log << std::endl;            
        }   

        log << std::endl;

        for (auto it : models)
        {
            delete it.second;
        }
	}
}