#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "../include/MeshConnectivity.h"
#include "testcases/HalfCylinder.h"
#include "testcases/Sphere.h"
#include "testcases/HalfTorus.h"
#include "testcases/Helix.h"
#include "testcases/Vee.h"
#include "igl/readOBJ.h"
#include <set>
#include <vector>
#include <map>
#include "igl/writePLY.h"
#include "igl/principal_curvature.h"

#include "models/ExactKichhoffLove.h"
#include "models/ExactMeanCurvatureSquared.h"
#include "models/S1Tan.h"
#include "models/S1Sin.h"
#include "models/S1Theta.h"
#include "models/SixVertexMidedgeAverage.h"
#include "models/VertexBasedQuadraticBending.h"
#include "models/IntrinsicVertexBasedQuadraticBending.h"
#include "models/CrouzeixRaviartQuadraticBending.h"
#include "models/DiscreteShells.h"
#include "models/DiscreteShellsTan.h"
#include "models/DKTBending.h"
#include "models/SecondOrderFEM.h"
#include "models/CorotationalHingeES.h"
#include "models/CorotationalHingeFS.h"

int main(int argc, char* argv[])
{
    double cylinderRadius = 0.0325;
    double cylinderHeight = 0.122;

    double veeWidth = 0.01;
    double veeHeight = 0.1;
	double PI = 3.1415926535898;
    std::vector<double> veeAngles = { 3.0*PI/4.0, PI/2.0,  3.0*PI/8.0, PI/4.0, 3.0*PI/16.0, PI/8.0, 3.0*PI/32.0, PI/16.0, PI/32.0 };

    double sphereRadius = 0.03;

    double outerTorusRadius = 0.03;
	double innerTorusRadius = 0.015;

	double helixRadius = 0.02;
	double helixHeight = 0.1;
	double helixSurfaceHeight = 0.1;

    double baseTriangleArea = 1e-8;

    // set up material parameters
    double thickness = 0.00010;

    double young = 1e11; // doesn't matter for static solves
    double poisson = 0.3;
    double lameAlpha = young * poisson / (1.0 - poisson * poisson);
    double lameBeta = young / 2.0 / (1.0 + poisson);

    std::map<std::string, Model*> models;
    models["Analytic (Koiter)"] = new ExactKirchhoffLove();
	models["Analytic (Squared H)"] = new ExactMeanCurvatureSquared();
    models["BAC"] = new S1Tan();
    models["DCS"] = new S1Sin();
    models["DSO"] = new S1Theta();
    models["Md.Avg."] = new SixVertexMidedgeAverage();
    models["QB(PL)"] = new VertexBasedQuadraticBending();
    models["QB(int)"] = new IntrinsicVertexBasedQuadraticBending();
	models["QB(CR)"] = new CrouzeixRaviartQuadraticBending();
	models["DS(theta)"] = new DiscreteShells();
    models["DS(tan)"] = new DiscreteShellsTan();
	models["DKT"] = new DKTBending();
    models["CHES"] = new CorotationalHingeES();
    models["CHFS"] = new CorotationalHingeFS();
    //models["DS(2nd)"] = new SecondOrderFEM();
    int steps = 5;
    double multiplier = 4;
    std::ofstream log("log.txt");
    log << std::setprecision(12);

	bool runCylinder = true;
    bool runVee = true;
	bool runSphere = true;
	bool runTorus = true;
    bool runHelix = true;

    bool genMeshes = true;


    if (runCylinder)
    {        
        for (bool regularMesh : {true, false})
        {
            log << "=========================" << std::endl;
            log << "Half-Cylinder, " << (regularMesh ? "regular" : "irregular") << " mesh" << std::endl;
            log << "=========================" << std::endl;
            log << "Vertices \\ ";
            bool first = true;
            for (auto model : models)
            {
                if (!first)
                    log << ", ";
                first = false;
                log << model.first;
            }
            log << std::endl;

            Eigen::MatrixXd origV;
            Eigen::MatrixXd rolledV;
            Eigen::MatrixXd origEdgeV;
            Eigen::MatrixXd rolledEdgeV;
            LibShell::MeshConnectivity mesh;
            double triangleArea = baseTriangleArea;
            
            makeHalfCylinder(regularMesh, false, cylinderRadius, cylinderHeight, 1.0, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);
            if (genMeshes)
            {
                igl::writePLY("halfcylinder-" + std::string(regularMesh ? "regular-" : "irregular-") + std::to_string(origV.rows()) + ".ply", rolledV, mesh.faces());
                igl::writePLY("halfcylinder-" + std::string(regularMesh ? "regular-flat" : "irregular-flat") + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());
            }

            for (int step = 0; step < steps; step++)
            {
                log << origV.rows() << ": ";
                bool first = true;
                for (auto model : models)
                {
                    if (!first)
                        log << ", ";
                    first = false;
                    log << model.second->measureHalfCylinderEnergy(mesh, origV, origEdgeV, rolledV, rolledEdgeV, thickness, lameAlpha, lameBeta, cylinderRadius, cylinderHeight);
                }
                log << std::endl;

                triangleArea *= multiplier;
                makeHalfCylinder(regularMesh, false, cylinderRadius, cylinderHeight, 1.0, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);                
                if (genMeshes)
                {
                    igl::writePLY("halfcylinder-" + std::string(regularMesh ? "regular-" : "irregular-") + std::to_string(origV.rows()) + ".ply", rolledV, mesh.faces());
                    igl::writePLY("halfcylinder-" + std::string(regularMesh ? "regular-flat" : "irregular-flat-") + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());
                }
            }
        }
        {
            log << "=========================" << std::endl;
            log << "Half-Cylinder, anisotropic mesh" << std::endl;
            log << "=========================" << std::endl;
            log << "Vertices \\ ";
            bool first = true;
            for (auto model : models)
            {
                if (!first)
                    log << ", ";
                first = false;
                log << model.first;
            }
            log << std::endl;

            Eigen::MatrixXd origV;
            Eigen::MatrixXd rolledV;
            Eigen::MatrixXd origEdgeV;
            Eigen::MatrixXd rolledEdgeV;
            LibShell::MeshConnectivity mesh;
            double triangleArea = baseTriangleArea;
            makeHalfCylinder(false, false, cylinderRadius, cylinderHeight, 4.0, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);
            if (genMeshes)
            {
                igl::writePLY("halfcylinder-aniso-" + std::to_string(origV.rows()) + ".ply", rolledV, mesh.faces());
                igl::writePLY("halfcylinder-aniso-flat-" + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());
            }
            for (int step = 0; step < steps; step++)
            {
                log << origV.rows() << ": ";
                bool first = true;
                for (auto model : models)
                {
                    if (!first)
                        log << ", ";
                    first = false;
                    log << model.second->measureHalfCylinderEnergy(mesh, origV, origEdgeV, rolledV, rolledEdgeV, thickness, lameAlpha, lameBeta, cylinderRadius, cylinderHeight);
                }
                log << std::endl;

                triangleArea *= multiplier;
                makeHalfCylinder(false, false, cylinderRadius, cylinderHeight, 4.0, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);
                if (genMeshes)
                {
                    igl::writePLY("halfcylinder-aniso-" + std::to_string(origV.rows()) + ".ply", rolledV, mesh.faces());
                    igl::writePLY("halfcylinder-aniso-flat-" + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());
                }
            }
        }
        
        {
            log << "=========================" << std::endl;
            log << "Rest-Curved Half-Cylinder, irregular mesh" << std::endl;
            log << "=========================" << std::endl;
            log << "Vertices \\ ";
            bool first = true;
            for (auto model : models)
            {
                if (!first)
                    log << ", ";
                first = false;
                log << model.first;
            }
            log << std::endl;

            Eigen::MatrixXd origV;
            Eigen::MatrixXd rolledV;
            Eigen::MatrixXd origEdgeV;
            Eigen::MatrixXd rolledEdgeV;
            LibShell::MeshConnectivity mesh;
            double triangleArea = baseTriangleArea;

            makeHalfCylinder(false, true, cylinderRadius, cylinderHeight, 1.0, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);
            if (genMeshes)
            {
                igl::writePLY("halfcylinder-restcurved-" + std::to_string(origV.rows()) + ".ply", rolledV, mesh.faces());
                igl::writePLY("halfcylinder-restcurved-flat-" + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());
            }

            for (int step = 0; step < steps; step++)
            {
                log << origV.rows() << ": ";
                bool first = true;
                for (auto model : models)
                {
                    if (!first)
                        log << ", ";
                    first = false;
                    log << model.second->measureRestCurvedHalfCylinderEnergy(mesh, origV, origEdgeV, rolledV, rolledEdgeV, thickness, lameAlpha, lameBeta, cylinderRadius, cylinderHeight);
                }
                log << std::endl;

                triangleArea *= multiplier;
                makeHalfCylinder(false, true, cylinderRadius, cylinderHeight, 1.0, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);
                if (genMeshes)
                {
                    igl::writePLY("halfcylinder-restcurved-" + std::to_string(origV.rows()) + ".ply", rolledV, mesh.faces());
                    igl::writePLY("halfcylinder-restcurved-flat-" + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());
                }
            }
        }
    }

    if (runVee)
    {
        int veeidx = 0;
        for (auto angle : veeAngles)
        {
            log << "=========================" << std::endl;
            log << "Vee, angle = " << angle << std::endl;
            log << "=========================" << std::endl;
            log << "Vertices \\ ";
            bool first = true;
            for (auto model : models)
            {
                if (!first)
                    log << ", ";
                first = false;
                log << model.first;
            }
            log << std::endl;

            Eigen::MatrixXd origV;
            Eigen::MatrixXd rolledV;
            Eigen::MatrixXd origEdgeV;
            Eigen::MatrixXd rolledEdgeV;
            LibShell::MeshConnectivity mesh;
            double triangleArea = baseTriangleArea;

            makeVee(veeWidth, veeHeight, angle, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);
            if (genMeshes)
            {
                std::stringstream ss;
                ss << "vee-" << veeidx << "-" << origV.rows() << ".ply";
                igl::writePLY(ss.str(), rolledV, mesh.faces());
                igl::writePLY("vee-flat-" + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());
            }

            for (int step = 0; step < steps; step++)
            {
                log << origV.rows() << ": ";
                bool first = true;
                for (auto model : models)
                {
                    if (!first)
                        log << ", ";
                    first = false;
                    log << model.second->measureVeeEnergy(mesh, origV, origEdgeV, rolledV, rolledEdgeV, thickness, lameAlpha, lameBeta, cylinderRadius, cylinderHeight);
                }
                log << std::endl;

                triangleArea *= multiplier;
                makeVee(veeWidth, veeHeight, angle, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);
                if (genMeshes)
                {
                    std::stringstream ss;
                    ss << "vee-" << veeidx << "-" << origV.rows() << ".ply";
                    igl::writePLY(ss.str(), rolledV, mesh.faces());
                    igl::writePLY("vee-flat-" + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());
                }
            }

            veeidx++;
        }
    }

    if (runSphere)
    {
        log << "=========================" << std::endl;
        log << "Sphere, irregular mesh" << std::endl;
        log << "=========================" << std::endl;
        log << "Vertices \\ ";
        bool first = true;
        for (auto model : models)
        {
            if (!first)
                log << ", ";
            first = false;
            log << model.first;
        }
        log << std::endl;

        double triangleArea = baseTriangleArea;
        Eigen::MatrixXd V;
        Eigen::MatrixXd edgeV;
        LibShell::MeshConnectivity mesh;

        makeSphere(sphereRadius, triangleArea, V, edgeV, mesh);        
        if (genMeshes)
        {
            igl::writePLY("sphere-" + std::to_string(V.rows()) + ".ply", V, mesh.faces());
        }

        for (int step = 0; step < steps; step++)
        {
            log << V.rows() << ": ";
            bool first = true;
            for (auto model : models)
            {
                if (!first)
                    log << ", ";
                first = false;
                log << model.second->measureSphereEnergy(mesh, V, edgeV, thickness, lameAlpha, lameBeta, sphereRadius);
            }
            log << std::endl;

            triangleArea *= multiplier;
            makeSphere(sphereRadius, triangleArea, V, edgeV, mesh);            
            if (genMeshes)
            {
                igl::writePLY("sphere-" + std::to_string(V.rows()) + ".ply", V, mesh.faces());
            }
        }
    }

    if (runTorus)
    {
        for (bool regularMesh : {true, false})
        {
            log << "=========================" << std::endl;
            log << "Half-Torus, " << (regularMesh ? "regular" : "irregular") << " mesh" << std::endl;
            log << "=========================" << std::endl;
            log << "Vertices \\ ";
            bool first = true;
            for (auto model : models)
            {
                if (!first)
                    log << ", ";
                first = false;
                log << model.first;
            }
            log << std::endl;

            double triangleArea = baseTriangleArea;
            Eigen::MatrixXd V;
            Eigen::MatrixXd edgeV;
            LibShell::MeshConnectivity mesh;

            makeHalfTorus(regularMesh, innerTorusRadius, outerTorusRadius, triangleArea, V, edgeV, mesh);
            if (genMeshes)
            {
                igl::writePLY("halftorus-" + std::string(regularMesh ? "regular-" : "irregular-") + std::to_string(V.rows()) + ".ply", V, mesh.faces());
            }

            for (int step = 0; step < steps; step++)
            {
                log << V.rows() << ": ";
                bool first = true;
                for (auto model : models)
                {
                    if (!first)
                        log << ", ";
                    first = false;
                    log << model.second->measureHalfTorusEnergy(mesh, V, edgeV, thickness, lameAlpha, lameBeta, innerTorusRadius, outerTorusRadius);
                }
                log << std::endl;

                triangleArea *= multiplier;
                makeHalfTorus(regularMesh, innerTorusRadius, outerTorusRadius, triangleArea, V, edgeV, mesh);                
                if (genMeshes)
                {
                    igl::writePLY("halftorus-" + std::string(regularMesh ? "regular-" : "irregular-") + std::to_string(V.rows()) + ".ply", V, mesh.faces());
                }
            }
        }
    }

    if (runHelix)
    {
        for (bool regularMesh : {true, false})
        {
            log << "=========================" << std::endl;
            log << "Helix Tangent Developable, " << (regularMesh ? "regular" : "irregular") << " mesh" << std::endl;
            log << "=========================" << std::endl;
            log << "Vertices \\ ";
            bool first = true;
            for (auto model : models)
            {
                if (!first)
                    log << ", ";
                first = false;
                log << model.first;
            }
            log << std::endl;

            Eigen::MatrixXd origV;
            Eigen::MatrixXd origEdgeV;
            Eigen::MatrixXd rolledV;
            Eigen::MatrixXd rolledEdgeV;
            LibShell::MeshConnectivity mesh;
            double triangleArea = baseTriangleArea;

            makeHelix(regularMesh, helixRadius, helixHeight, helixSurfaceHeight, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);
            if (genMeshes)
            {
                igl::writePLY("helix-" + std::string(regularMesh ? "regular-" : "irregular-") + std::to_string(origV.rows()) + ".ply", rolledV, mesh.faces());
                igl::writePLY("helix-" + std::string(regularMesh ? "regular-flat-" : "irregular-flat-") + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());

            }

            for (int step = 0; step < steps; step++)
            {
                log << origV.rows() << ": ";
                bool first = true;
                for (auto model : models)
                {
                    if (!first)
                        log << ", ";
                    first = false;
					log << model.second->measureHelixEnergy(mesh, origV, origEdgeV, rolledV, rolledEdgeV, thickness, lameAlpha, lameBeta, helixRadius, helixHeight, helixSurfaceHeight);
                }
                log << std::endl;

                triangleArea *= multiplier;
                makeHelix(regularMesh, helixRadius, helixHeight, helixSurfaceHeight, triangleArea, origV, origEdgeV, rolledV, rolledEdgeV, mesh);                
                if (genMeshes)
                {
                    igl::writePLY("helix-" + std::string(regularMesh ? "regular-" : "irregular-") + std::to_string(origV.rows()) + ".ply", rolledV, mesh.faces());
                    igl::writePLY("helix-" + std::string(regularMesh ? "regular-flat-" : "irregular-flat-") + std::to_string(origV.rows()) + ".ply", origV, mesh.faces());

                }
            }
        }
    }
}
