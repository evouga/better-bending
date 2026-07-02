#include "Helix.h"
#include "igl/triangle/triangulate.h"
#include <sstream>
#include <cassert>
#include <iomanip>
#include "igl/remove_unreferenced.h"
#include "igl/boundary_loop.h"
#include "../include/MeshConnectivity.h"

void makeHelix(bool regular, double helixRadius, double helixHeight, double surfaceHeight, double triangleArea,
    Eigen::MatrixXd& flatV,
    Eigen::MatrixXd& flatEdgeV,
    Eigen::MatrixXd& V,
    Eigen::MatrixXd& edgeV,
    LibShell::MeshConnectivity &mesh)
{
    constexpr double PI = 3.1415926535898;
    double targetlength = std::sqrt(2.0 * triangleArea);
    
	double alpha = std::sqrt(helixRadius * helixRadius + helixHeight * helixHeight / 4.0 / PI / PI);
	
    int W = std::max(1, int(2.0 * PI * helixRadius / targetlength));
    int H = std::max(1, int(surfaceHeight / targetlength));

    if (regular)
    {
        flatV.resize((W + 1) * (H + 1), 3);
        V.resize((W + 1) * (H + 1), 3);
        Eigen::MatrixXi F(2 * W * H, 3);
        Eigen::MatrixXd paramV((W + 1) * (H + 1), 2);
        int curface = 0;
        for (int i = 0; i <= H; i++)
        {
            for (int j = 0; j <= W; j++)
            {
                int idx = i * (W + 1) + j;
                double x = 2.0 * PI * double(j) / double(W);
                double y = surfaceHeight / helixRadius * (1.0 + double(i) / double(H) );
				flatV(idx, 0) = alpha * alpha / helixRadius * std::cos(helixRadius / alpha * x) - y * alpha * std::sin(helixRadius / alpha * x);
				flatV(idx, 1) = alpha * alpha / helixRadius * std::sin(helixRadius / alpha * x) + y * alpha * std::cos(helixRadius / alpha * x);
                flatV(idx, 2) = 0;
                V(idx, 0) = helixRadius * std::cos(x) - y * helixRadius * std::sin(x);
				V(idx, 1) = helixRadius * std::sin(x) + y * helixRadius * std::cos(x);
				V(idx, 2) = helixHeight * (x + y) / (2.0 * PI);
                paramV(idx, 0) = x;
                paramV(idx, 1) = y;
                if (i > 0 && j > 0)
                {
                    int idxm1m1 = (i - 1) * (W + 1) + (j - 1);
                    int idxm1m0 = (i - 1) * (W + 1) + j;
                    F(curface, 0) = idxm1m1;
                    F(curface, 1) = idxm1m0;
                    F(curface, 2) = idx;
                    int idxm0m1 = i * (W + 1) + (j - 1);
                    F(curface + 1, 0) = idxm1m1;
                    F(curface + 1, 1) = idx;
                    F(curface + 1, 2) = idxm0m1;
                    curface += 2;
                }
            }
        }
        mesh = LibShell::MeshConnectivity(F);
        int nedges = mesh.nEdges();
        flatEdgeV.resize(nedges, 3);
        edgeV.resize(nedges, 3);
        for (int i = 0; i < nedges; i++)
        {
            int v0 = mesh.edgeVertex(i, 0);
            int v1 = mesh.edgeVertex(i, 1);
            Eigen::Vector2d midpt = 0.5 * (paramV.row(v0) + paramV.row(v1)).transpose();
            double x = midpt[0];
            double y = midpt[1];
            flatEdgeV(i, 0) = alpha * alpha / helixRadius * std::cos(helixRadius / alpha * x) - y * alpha * std::sin(helixRadius / alpha * x);
            flatEdgeV(i, 1) = alpha * alpha / helixRadius * std::sin(helixRadius / alpha * x) + y * alpha * std::cos(helixRadius / alpha * x);
            flatEdgeV(i, 2) = 0;
            edgeV(i, 0) = helixRadius * std::cos(x) - y * helixRadius * std::sin(x);
            edgeV(i, 1) = helixRadius * std::sin(x) + y * helixRadius * std::cos(x);
            edgeV(i, 2) = helixHeight * (x + y) / (2.0 * PI);
        }
    }
    else
    {
        Eigen::MatrixXd Vin(2 * W + 2 * H, 2);
        Eigen::MatrixXi E(2 * W + 2 * H, 2);
        Eigen::MatrixXd dummyH(0, 2);
        Eigen::MatrixXd V2;
        Eigen::MatrixXi F2;
        
        int vrow = 0;
        int erow = 0;
        // top boundary
        for (int i = 1; i < W; i++)
        {
            Vin(vrow, 0) = double(i) / double(W) * 2.0 * PI * helixRadius;
            Vin(vrow, 1) = surfaceHeight;
            if (i > 1)
            {
                E(erow, 0) = vrow - 1;
                E(erow, 1) = vrow;
                erow++;
            }
            vrow++;
        }
        // bottom boundary
        for (int i = 1; i < W; i++)
        {
            Vin(vrow, 0) = double(i) / double(W) * 2.0 * PI * helixRadius;
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
        for (int i = 0; i <= H; i++)
        {
            Vin(vrow, 0) = 0;
            Vin(vrow, 1) = double(i) / double(H) * surfaceHeight;
            if (i > 0)
            {
                E(erow, 0) = vrow - 1;
                E(erow, 1) = vrow;
                erow++;
            }
            vrow++;
        }
        // right boundary    
        for (int i = 0; i <= H; i++)
        {
            Vin(vrow, 0) = 2.0 * PI * helixRadius;
            Vin(vrow, 1) = double(i) / double(H) * surfaceHeight;
            if (i > 0)
            {
                E(erow, 0) = vrow - 1;
                E(erow, 1) = vrow;
                erow++;
            }
            vrow++;
        }
        // missing four edges
        E(erow, 0) = (W - 1) - 1;
        E(erow, 1) = 2 * (W - 1) + 2 * (H + 1) - 1;
        erow++;
        E(erow, 0) = 2 * (W - 1) + (H + 1);
        E(erow, 1) = 2 * (W - 1) - 1;
        erow++;
        E(erow, 0) = W - 1;
        E(erow, 1) = 2 * (W - 1);
        erow++;
        E(erow, 0) = 2 * (W - 1) + (H + 1) - 1;
        E(erow, 1) = 0;
        erow++;

        assert(vrow == 2 * H + 2 * W);
        assert(erow == 2 * H + 2 * W);
        std::stringstream ss;
        ss << "a" << std::setprecision(30) << std::fixed << triangleArea << "qDY";
        igl::triangle::triangulate(Vin, E, dummyH, ss.str(), V2, F2);

        // roll up

        int nverts = V2.rows();

        flatV.resize(nverts, 3);
        Eigen::MatrixXd rolledV(nverts, 3);
        Eigen::MatrixXd paramV(nverts, 2);

        for (int i = 0; i < nverts; i++)
        {
            Eigen::Vector2d q = V2.row(i).transpose();
            double x = q[0] / helixRadius;
			double y = (surfaceHeight + q[1]) / helixRadius;
            flatV(i, 0) = alpha * alpha / helixRadius * std::cos(helixRadius / alpha * x) - y * alpha * std::sin(helixRadius / alpha * x);
            flatV(i, 1) = alpha * alpha / helixRadius * std::sin(helixRadius / alpha * x) + y * alpha * std::cos(helixRadius / alpha * x);
            flatV(i, 2) = 0;
            rolledV(i, 0) = helixRadius * std::cos(x) - y * helixRadius * std::sin(x);
            rolledV(i, 1) = helixRadius * std::sin(x) + y * helixRadius * std::cos(x);
            rolledV(i, 2) = helixHeight * (x + y) / (2.0 * PI);
            paramV(i, 0) = x;
            paramV(i, 1) = y;
        }

        V = rolledV;
        mesh = LibShell::MeshConnectivity(F2);
        int nedges = mesh.nEdges();
        flatEdgeV.resize(nedges, 3);
        edgeV.resize(nedges, 3);
        for (int i = 0; i < nedges; i++)
        {
            int v0 = mesh.edgeVertex(i, 0);
            int v1 = mesh.edgeVertex(i, 1);
            Eigen::Vector2d midpt = 0.5 * (paramV.row(v0) + paramV.row(v1)).transpose();
            double x = midpt[0];
            double y = midpt[1];
            flatEdgeV(i, 0) = alpha * alpha / helixRadius * std::cos(helixRadius / alpha * x) - y * alpha * std::sin(helixRadius / alpha * x);
            flatEdgeV(i, 1) = alpha * alpha / helixRadius * std::sin(helixRadius / alpha * x) + y * alpha * std::cos(helixRadius / alpha * x);
            flatEdgeV(i, 2) = 0;
            edgeV(i, 0) = helixRadius * std::cos(x) - y * helixRadius * std::sin(x);
            edgeV(i, 1) = helixRadius * std::sin(x) + y * helixRadius * std::cos(x);
            edgeV(i, 2) = helixHeight * (x + y) / (2.0 * PI);
        }
    }
}
