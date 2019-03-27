/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include <minmax.h>

NORI_NAMESPACE_BEGIN


Point2f Warp::squareToUniformSquare(const Point2f &sample)
{
	return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample)
{
	return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

float nori::Warp::inverseCumulateTent(const float& v)
{
	return v <= 0.5f ? sqrt(2 * v) - 1 : 1 - sqrt(2 * (1 - v));
}

float nori::Warp::tent(const float& v)
{
	return (v >= -1 && v <= 1) ? 1 - abs(v) : 0;
}

Point2f Warp::squareToTent(const Point2f &sample)
{
	return Point2f(inverseCumulateTent(sample[0]), inverseCumulateTent(sample[1]));
}

float Warp::squareToTentPdf(const Point2f &p)
{
	return tent(p[0]) * tent(p[1]);
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) 
{
	float r = sqrt(sample.x());
	float theta = 2.0f * M_PI * sample.y();

	return Point2f(r * cos(theta), r * sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p)
{
	return (p.x() * p.x() + p.y() * p.y() <= 1.0f) ? 1 / M_PI : 0.0f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample)
{
	float z = 1 - 2 * sample.x();

	float theta = 2 * M_PI * sample.y();
	float phi = 1 - z * z;

	float x = cos(theta) * sqrt(phi);
	float y = sin(theta) * sqrt(phi);

	return Vector3f(x, y, z);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v)
{
	return INV_FOURPI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample)
{
	float theta = 2 * M_PI * sample.y();
	float phi = 1 - sample.x() * sample.x();

	float z = sample.x();
	float x = cos(theta) * sqrt(phi);
	float y = sin(theta) * sqrt(phi);

	return Vector3f(x, y, z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v)
{
	if (v.z() < 0)
		return 0;
	return INV_TWOPI;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample)
{
	Point2f disk = squareToUniformDisk(sample);

	Vector3f ret = Vector3f(disk.x(), disk.y(), 0);
	ret[2] = sqrt(max(0.0f, 1.0f - ret.x() * ret.x() - ret.y() * ret.y()));
	return ret;
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v)
{
	if (v.z() < 0)
		return 0;
	return v.z() * INV_PI;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha)
{
	float phi = 2 * M_PI * sample.x();
	float theta = atan(sqrt(-alpha * alpha * log(1 - sample.y())));

	float x = sin(theta) * cos(phi);
	float y = sin(theta) * sin(phi);
	float z = cos(theta);

	return Vector3f(x, y, z);
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha)
{
	if (m.z() <= 0)
		return 0;

	float alpha_2 = alpha * alpha;
	float tan_theta_2 = (m.x() * m.x() + m.y() * m.y()) / (m.z() * m.z());
	float cos_theta_3 = m.z() * m.z() * m.z();

	float azimuthal = INV_TWOPI;
	float longitudinal = (2 * exp(-tan_theta_2 / alpha_2)) / (alpha_2 * cos_theta_3);

	return azimuthal * longitudinal;
}


NORI_NAMESPACE_END
