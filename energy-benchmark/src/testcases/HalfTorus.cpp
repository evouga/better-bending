#include "HalfTorus.h"
#include "igl/triangle/triangulate.h"
#include <sstream>
#include <cassert>
#include <iomanip>
#include "igl/remove_unreferenced.h"
#include "igl/boundary_loop.h"
#include "../include/MeshConnectivity.h"

void makeHalfTorus(bool regular, double innerRadius, double outerRadius, double triangleArea,
    Eigen::MatrixXd& V,
    Eigen::MatrixXd &edgeV,
    LibShell::MeshConnectivity &mesh)
{
    constexpr double PI = 3.1415926535898;
    double targetlength = std::sqrt(2.0 * triangleArea);

    int W = std::max(1, int(2.0 * PI * outerRadius / targetlength));
    int H = std::max(1, int(PI * innerRadius / targetlength));

    if (regular)
    {
        V.resize(W * (H + 1), 3);
        Eigen::MatrixXd paramV(W * (H + 1), 2);
        Eigen::MatrixXi F(2 * W * H, 3);
        int curface = 0;
        for (int i = 0; i <= H; i++)
        {
            for (int j = 0; j < W; j++)
            {
                int idx = i * W + j;
                double x = double(j) / double(W);
                double y = double(i) / double(H);
                paramV(idx, 0) = x;
                paramV(idx, 1) = y;
				V(idx, 0) = (outerRadius + innerRadius * std::cos(y * PI)) * std::cos(x * 2.0 * PI);
				V(idx, 1) = (outerRadius + innerRadius * std::cos(y * PI)) * std::sin(x * 2.0 * PI);
				V(idx, 2) = innerRadius * std::sin(y * PI);
                if (i > 0 && j > 0)
                {
                    int idxm1m1 = (i - 1) * W + (j - 1);
                    int idxm1m0 = (i - 1) * W + j;
                    F(curface, 0) = idxm1m1;
                    F(curface, 1) = idxm1m0;
                    F(curface, 2) = idx;
                    int idxm0m1 = i * W + (j - 1);
                    F(curface + 1, 0) = idxm1m1;
                    F(curface + 1, 1) = idx;
                    F(curface + 1, 2) = idxm0m1;
                    curface += 2;
                }
            }
        }
        for (int i = 0; i < H; i++)
        {
            int idxm1 = i * W;
            int idxm0 = (i + 1) * W;
            F(curface, 0) = idxm0 + W - 1;
            F(curface, 1) = idxm1 + W - 1;
            F(curface, 2) = idxm0;
            curface++;
            F(curface, 0) = idxm0;
            F(curface, 1) = idxm1 + W - 1;
            F(curface, 2) = idxm1;
			curface++;
        }
        mesh = LibShell::MeshConnectivity(F);
        int nedges = mesh.nEdges();
        edgeV.resize(nedges, 3);

        for (int i = 0; i < nedges; i++)
        {
            int v0 = mesh.edgeVertex(i, 0);
            int v1 = mesh.edgeVertex(i, 1);
			Eigen::Vector2d pt0 = paramV.row(v0).transpose();
			Eigen::Vector2d pt1 = paramV.row(v1).transpose();
            if(std::fabs(pt0[0]-pt1[0]) > 0.5)
            {
                // edge crosses the seam
                if(pt0[0] > pt1[0])
                {
                    std::swap(pt0, pt1);
                }
                pt0[0] += 1.0;                
			}
            Eigen::Vector2d midpt = 0.5 * (pt0 + pt1);
            double x = midpt[0];
            double y = midpt[1];
            edgeV(i, 0) = (outerRadius + innerRadius * std::cos(y * PI)) * std::cos(x * 2.0 * PI);
            edgeV(i, 1) = (outerRadius + innerRadius * std::cos(y * PI)) * std::sin(x * 2.0 * PI);
            edgeV(i, 2) = innerRadius * std::sin(y * PI);
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
            Vin(vrow, 0) = double(i) / double(W) * 2.0 * PI * outerRadius;
            Vin(vrow, 1) = PI * innerRadius;
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
            Vin(vrow, 0) = double(i) / double(W) * 2.0 * PI * outerRadius;
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
            Vin(vrow, 1) = double(i) / double(H) * PI * innerRadius;
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
            Vin(vrow, 0) = 2.0 * PI * outerRadius;
            Vin(vrow, 1) = double(i) / double(H) * PI * innerRadius;
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

        Eigen::MatrixXd rolledV(nverts, 3);
        Eigen::MatrixXd paramV(nverts, 2);

        for (int i = 0; i < nverts; i++)
        {
            Eigen::Vector2d q = V2.row(i).transpose();
            Eigen::Vector3d rolledq;
			double x = q[0] / (2.0 * PI * outerRadius);
			double y = q[1] / (PI * innerRadius);
            paramV(i, 0) = x;
            paramV(i, 1) = y;
            rolledq[0] = (outerRadius + innerRadius * std::cos(y * PI)) * std::cos(x * 2.0 * PI);
            rolledq[1] = (outerRadius + innerRadius * std::cos(y * PI)) * std::sin(x * 2.0 * PI);
            rolledq[2] = innerRadius * std::sin(y * PI);
            rolledV.row(i) = rolledq.transpose();
        }

        std::map<int, int> vertexMap;
        for (int i = 0; i <= H; i++)
        {
            vertexMap[2 * (W - 1) + i + H + 1] = 2 * (W - 1) + i;
        }
        for (int i = 0; i < F2.rows(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
				auto it = vertexMap.find(F2(i, j));
                if(it != vertexMap.end())
                {
                    F2(i, j) = it->second;
				}
            }
        }
        Eigen::VectorXd I;
        Eigen::MatrixXi F;
        Eigen::MatrixXd prunedV;
        igl::remove_unreferenced(rolledV, F2, V, F, I);
        igl::remove_unreferenced(paramV, F2, prunedV, F, I);
        mesh = LibShell::MeshConnectivity(F);
        int nedges = mesh.nEdges();
        edgeV.resize(nedges, 3);

        for (int i = 0; i < nedges; i++)
        {
            int v0 = mesh.edgeVertex(i, 0);
            int v1 = mesh.edgeVertex(i, 1);
            Eigen::Vector2d pt0 = prunedV.row(v0).transpose();
            Eigen::Vector2d pt1 = prunedV.row(v1).transpose();
            if (std::fabs(pt0[0] - pt1[0]) > 0.5)
            {
                // edge crosses the seam
                if (pt0[0] > pt1[0])
                {
                    std::swap(pt0, pt1);
                }
                pt0[0] += 1.0;
            }
            Eigen::Vector2d midpt = 0.5 * (pt0 + pt1);
            double x = midpt[0];
            double y = midpt[1];
            edgeV(i, 0) = (outerRadius + innerRadius * std::cos(y * PI)) * std::cos(x * 2.0 * PI);
            edgeV(i, 1) = (outerRadius + innerRadius * std::cos(y * PI)) * std::sin(x * 2.0 * PI);
            edgeV(i, 2) = innerRadius * std::sin(y * PI);
        }
    }
}