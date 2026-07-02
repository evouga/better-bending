#include "Vee.h"
#include "igl/triangle/triangulate.h"
#include <sstream>
#include <cassert>
#include <iomanip>
#include "igl/remove_unreferenced.h"
#include "igl/boundary_loop.h"
#include "../include/MeshConnectivity.h"

void makeVee(double width, double height, double creaseAngle, double triangleArea,
    Eigen::MatrixXd& flatV,
    Eigen::MatrixXd& flatEdgeV,
    Eigen::MatrixXd& V,
    Eigen::MatrixXd& edgeV,
    LibShell::MeshConnectivity& mesh)
{
    double targetlength = std::sqrt(2.0 * triangleArea);

    int W = 2 * std::max(1, int(width / targetlength / 2.0));
    int H = std::max(1, int(height / targetlength));

	double creaseX = 0.5 * width;

    flatV.resize((W + 1) * (H + 1), 3);
    V.resize((W + 1) * (H + 1), 3);
    Eigen::MatrixXi F(2 * W * H, 3);
    int curface = 0;
    for (int i = 0; i <= H; i++)
    {
        for (int j = 0; j <= W; j++)
        {
            int idx = i * (W + 1) + j;
            flatV(idx, 0) = double(j) / double(W) * width;
            flatV(idx, 1) = double(i) / double(H) * height;
            flatV(idx, 2) = 0;
			double disp = flatV(idx,0) - creaseX;
            V(idx, 0) = creaseX + disp * std::sin(creaseAngle / 2.0);
            V(idx, 1) = flatV(idx, 1);
            V(idx, 2) = std::fabs(disp) * std::cos(creaseAngle / 2.0);
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
        Eigen::Vector3d midpt = 0.5 * (flatV.row(v0) + flatV.row(v1)).transpose();
        double disp = midpt[0] - creaseX;
        flatEdgeV.row(i) = midpt.transpose();
		edgeV(i, 0) = creaseX + disp * std::sin(creaseAngle / 2.0);
        edgeV(i, 1) = midpt[1];
		edgeV(i, 2) = std::fabs(disp) * std::cos(creaseAngle / 2.0);
    }
}