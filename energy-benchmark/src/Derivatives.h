#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include <Eigen/Core>
#include <vector>

Eigen::Matrix3d crossMatrix(const Eigen::Vector3d& v)
{
	Eigen::Matrix3d m;
	m << 0.0, -v.z(), v.y(),
		v.z(), 0.0, -v.x(),
		-v.y(), v.x(), 0.0;
	return m;
}

Eigen::Vector3d weightedNormal(const Eigen::Vector3d& u, const Eigen::Vector3d& v,
	Eigen::Matrix<double, 3, 6>* deriv,
	std::vector<Eigen::Matrix<double, 6, 6> >* hess)
{
	if (deriv)
	{
		Eigen::Matrix3d cu = crossMatrix(u);
		Eigen::Matrix3d cv = crossMatrix(v);
		deriv->block<3, 3>(0, 0) = -cv;
		deriv->block<3, 3>(0, 3) = cu;
	}
	if (hess)
	{
		hess->resize(3);
		for (int i = 0; i < 3; i++)
		{
			(*hess)[i].setZero();
		}
		(*hess)[0](1, 5) = 1.0;
		(*hess)[0](2, 4) = -1.0;
		(*hess)[0](4, 2) = -1.0;
		(*hess)[0](5, 1) = 1.0;
		(*hess)[1](0, 5) = -1.0;
		(*hess)[1](2, 3) = 1.0;
		(*hess)[1](3, 2) = 1.0;
		(*hess)[1](5, 0) = -1.0;
		(*hess)[2](0, 4) = 1.0;
		(*hess)[2](1, 3) = -1.0;
		(*hess)[2](3, 1) = -1.0;
		(*hess)[2](4, 0) = 1.0;
	}
	return u.cross(v);
}

Eigen::Vector3d normal(const Eigen::Vector3d& u, const Eigen::Vector3d& v,
	Eigen::Matrix<double, 3, 6>* deriv,
	std::vector<Eigen::Matrix<double, 6, 6> >* hess)
{
	Eigen::Matrix<double, 3, 6> wderiv;
	std::vector<Eigen::Matrix<double, 6, 6> > whess;
	Eigen::Vector3d w = weightedNormal(u, v, deriv || hess ? &wderiv : NULL, hess ? &whess : NULL);
	double wnorm = w.norm();
	Eigen::Vector3d what = w / wnorm;
	if (deriv)
	{
		*deriv = 1.0 / wnorm * (Eigen::Matrix3d::Identity() - what * what.transpose()) * wderiv;
	}
	if (hess)
	{
		hess->resize(3);
		Eigen::Matrix<double, 6, 6> whathess;
		whathess.setZero();
		for(int i=0; i<3; i++)
			whathess += what[i] * whess[i];

		for (int i = 0; i < 3; i++)
		{
			Eigen::Vector3d ei = Eigen::Vector3d::Zero();
			ei[i] = 1.0;
			Eigen::Matrix3d term1mat = 1.0 / wnorm / wnorm * ((3.0 * what * what.transpose() - Eigen::Matrix3d::Identity()) * what[i] - what * ei.transpose() - ei * what.transpose());
			(*hess)[i] = wderiv.transpose() * term1mat * wderiv;
			(*hess)[i] += 1.0 / wnorm * (whess[i] - whathess * what[i]);
		}
	}
	return what;
}

template<int N> 
Eigen::Vector3d faceNormal(
	const Eigen::Matrix<double, 6, 1>& dr,
	const Eigen::Matrix<double, 6, N>& dxdr,
	Eigen::Matrix<double, 3, N>* deriv,
	std::vector<Eigen::Matrix<double, N, N> >* hess)
{
	Eigen::Matrix<double, 3, 6> nderiv;
	std::vector<Eigen::Matrix<double, 6, 6> > nhess;
	Eigen::Vector3d nhat = normal(dr.block<3, 1>(0, 0).transpose(), dr.block<3, 1>(3, 0).transpose(), (deriv || hess) ? &nderiv : NULL, hess ? &nhess : NULL);
	if (deriv)
	{
		*deriv = nderiv * dxdr;		
	}
	if (hess)
	{
		hess->resize(3);
		for (int i = 0; i < 3; i++)
		{
			(*hess)[i] = dxdr.transpose() * nhess[i] * dxdr;
		}
	}
	return nhat;
}

template<int N>
Eigen::Vector3d faceAreaNormal(
	const Eigen::Matrix<double, 6, 1>& dr,
	const Eigen::Matrix<double, 6, N>& dxdr,
	Eigen::Matrix<double, 3, N>* deriv,
	std::vector<Eigen::Matrix<double, N, N> >* hess)
{
	Eigen::Matrix<double, 3, 6> nderiv;
	std::vector<Eigen::Matrix<double, 6, 6> > nhess;
	Eigen::Vector3d nhat = weightedNormal(dr.block<3, 1>(0, 0).transpose(), dr.block<3, 1>(3, 0).transpose(), (deriv || hess) ? &nderiv : NULL, hess ? &nhess : NULL);
	if (deriv)
	{
		*deriv = nderiv * dxdr;
	}
	if (hess)
	{
		hess->resize(3);
		for (int i = 0; i < 3; i++)
		{
			(*hess)[i] = dxdr.transpose() * nhess[i] * dxdr;
		}
	}
	return nhat;
}

template<int N> 
Eigen::Vector4d secondFundamentalForm(
	const Eigen::Matrix<double, 6, 1>& dr, const std::vector<Eigen::Vector3d>& d2r,
	const Eigen::Matrix<double, 6, N>& dxdr, const std::vector < Eigen::Matrix<double, 3, N> >& dxd2r,
	Eigen::Matrix<double, 4, N>* deriv,
	std::vector<Eigen::Matrix<double, N, N> >* hess)
{
	Eigen::Matrix<double, 3, 6> nderiv;
	std::vector<Eigen::Matrix<double, 6, 6> > nhess;
	Eigen::Vector3d nhat = normal(dr.block<3, 1>(0, 0), dr.block<3, 1>(3, 0), (deriv || hess) ? &nderiv : NULL, hess ? &nhess : NULL);
	Eigen::Vector4d result;
	for (int i = 0; i < 4; i++)
	{
		result[i] = -d2r[i].dot(nhat);
	}
	if (deriv)
	{
		for (int i = 0; i < 4; i++)
		{
			deriv->row(i) = -nhat.transpose() * dxd2r[i] - d2r[i].transpose() * nderiv * dxdr;
		}
	}
	if (hess)
	{
		hess->resize(4);
		
		for (int i = 0; i < 4; i++)
		{
			(*hess)[i] = -dxdr.transpose() * nderiv.transpose() * dxd2r[i];
			(*hess)[i] += -dxd2r[i].transpose() * nderiv * dxdr;
	
			Eigen::Matrix<double, 6, 6> rihess;
			rihess.setZero();
			for (int j = 0; j < 3; j++)
			{
				rihess += d2r[i][j] * nhess[j];
			}

			(*hess)[i] += -dxdr.transpose() * rihess * dxdr;
		}
	}

	return result;
}

template<int N>
Eigen::Vector4d strainTensor(
	const Eigen::Matrix<double, 6, 1>& dr, const std::vector<Eigen::Vector3d>& d2r,
	const Eigen::Matrix<double, 6, N>& dxdr, const std::vector < Eigen::Matrix<double, 3, N> >& dxd2r,
	const Eigen::Vector4d& abar,
	const Eigen::Vector4d& bbar,
	Eigen::Matrix<double, 4, N>* deriv,
	std::vector<Eigen::Matrix<double, N, N> >* hess
	)
{
	Eigen::Matrix<double, 4, N> bderiv;
	std::vector<Eigen::Matrix<double, N, N> > bhess;
	Eigen::Vector4d b = secondFundamentalForm<N>(dr, d2r, dxdr, dxd2r, (deriv || hess) ? &bderiv : NULL, hess ? &bhess : NULL);
	double abardet = abar[0] * abar[3] - abar[1] * abar[2];
	Eigen::Vector4d bdiff = b - bbar;
	Eigen::Vector4d result;
	result << (abar[3] * bdiff[0] - abar[1] * bdiff[2]) / abardet,
		(abar[3] * bdiff[1] - abar[1] * bdiff[3]) / abardet,
		(-abar[2] * bdiff[0] + abar[0] * bdiff[2]) / abardet,
		(-abar[2] * bdiff[1] + abar[0] * bdiff[3]) / abardet;

	if (deriv)
	{
		deriv->row(0) = (abar[3] * bderiv.row(0) - abar[1] * bderiv.row(2)) / abardet;
		deriv->row(1) = (abar[3] * bderiv.row(1) - abar[1] * bderiv.row(3)) / abardet;
		deriv->row(2) = (-abar[2] * bderiv.row(0) + abar[0] * bderiv.row(2)) / abardet;
		deriv->row(3) = (-abar[2] * bderiv.row(1) + abar[0] * bderiv.row(3)) / abardet;
	}
	if (hess)
	{
		hess->resize(4);
		(*hess)[0] = (abar[3] * bhess[0] - abar[1] * bhess[2]) / abardet;
		(*hess)[1] = (abar[3] * bhess[1] - abar[1] * bhess[3]) / abardet;
		(*hess)[2] = (-abar[2] * bhess[0] + abar[0] * bhess[2]) / abardet;
		(*hess)[3] = (-abar[2] * bhess[1] + abar[0] * bhess[3]) / abardet;
	}
	return result;
}

template<int N>
double StVKEnergy(
	const Eigen::Matrix<double, 6, 1>& dr, const std::vector<Eigen::Vector3d>& d2r,
	const Eigen::Matrix<double, 6, N>& dxdr, const std::vector < Eigen::Matrix<double, 3, N> >& dxd2r,
	const Eigen::Vector4d& abar,
	const Eigen::Vector4d& bbar,
	double lameAlpha, double lameBeta,
	Eigen::Matrix<double, 1, N>* deriv,
	Eigen::Matrix<double, N, N>* hess)
{
	Eigen::Matrix<double, 4, N> strainDeriv;
	std::vector<Eigen::Matrix<double, N, N> > strainHess;
	Eigen::Vector4d strain = strainTensor<N>(dr, d2r, dxdr, dxd2r, abar, bbar,
		(deriv || hess) ? &strainDeriv : NULL,
		hess ? &strainHess : NULL);
	
	double energy = 0.5 * lameAlpha * (strain[0] + strain[3]) * (strain[0] + strain[3])
		+ lameBeta * (strain[0] * strain[0] + 2.0 * strain[2] * strain[1] + strain[3] * strain[3]);

	if (deriv)
	{
		(*deriv) = lameAlpha * (strain[0] + strain[3]) * (strainDeriv.row(0) + strainDeriv.row(3))
			+ 2.0 * lameBeta * (strain[0] * strainDeriv.row(0) + strain[2] * strainDeriv.row(1)
				+ strain[1] * strainDeriv.row(2) + strain[3] * strainDeriv.row(3));
	}

	if (hess)
	{
		(*hess) = lameAlpha * (strain[0] + strain[3]) * (strainHess[0] + strainHess[3])
			+ lameAlpha * (strainDeriv.row(0) + strainDeriv.row(3)).transpose() * (strainDeriv.row(0) + strainDeriv.row(3));
		(*hess) += 2.0 * lameBeta * (strain[0] * strainHess[0] + strain[2] * strainHess[1]
				+ strain[1] * strainHess[2] + strain[3] * strainHess[3]);
		(*hess) += 2.0 * lameBeta * (strainDeriv.row(0).transpose() * strainDeriv.row(0)
			+ strainDeriv.row(1).transpose() * strainDeriv.row(2)
			+ strainDeriv.row(2).transpose() * strainDeriv.row(1)
			+ strainDeriv.row(3).transpose() * strainDeriv.row(3));
	}

	return energy;
}

template<int N>
double MeanCurvatureSquaredEnergy(
	const Eigen::Matrix<double, 6, 1>& dr, const std::vector<Eigen::Vector3d>& d2r,
	const Eigen::Matrix<double, 6, N>& dxdr, const std::vector < Eigen::Matrix<double, 3, N> >& dxd2r,
	const Eigen::Vector4d& abar,
	const Eigen::Vector4d& bbar,
	double lameAlpha, double lameBeta,
	Eigen::Matrix<double, 1, N>* deriv,
	Eigen::Matrix<double, N, N>* hess
	)
{
	Eigen::Matrix<double, 4, N> strainDeriv;
	std::vector<Eigen::Matrix<double, N, N> > strainHess;
	Eigen::Vector4d strain = strainTensor<N>(dr, d2r, dxdr, dxd2r, abar, bbar,
		(deriv || hess) ? &strainDeriv : NULL,
		hess ? &strainHess : NULL
		);

	double energy = (0.5 * lameAlpha + lameBeta) * (strain[0] + strain[3]) * (strain[0] + strain[3]);

	if (deriv)
	{
		(*deriv) = (lameAlpha + 2.0 * lameBeta) * (strain[0] + strain[3]) * (strainDeriv.row(0) + strainDeriv.row(3));
	}

	if (hess)
	{
		Eigen::Matrix<double, 1, N> derivSum = strainDeriv.row(0) + strainDeriv.row(3);
		(*hess) = (lameAlpha + 2.0 * lameBeta) * (strain[0] + strain[3]) * (strainHess[0] + strainHess[3])
			+ (lameAlpha + 2.0 * lameBeta) * (derivSum.transpose() * derivSum);
	}
	return energy;
}

/*
 * Enforces the inequality constraint x >= 0, with zero force when x >= 1.
 */
double barrierPotential(double x, double& deriv, double& hess)
{
	deriv = 0;
	hess = 0;
	if (x <= 0)
		return std::numeric_limits<double>::infinity();
	double result = -0.5 * (x - 1.0) * (x - 3.0) - std::log(x);
	if(result < 0)
		return 0;
	deriv = -(x - 1.0) * (x - 1.0) / x;
	hess = (1.0 - x) * (1.0 + x) / x / x;
	return result;
}

#endif