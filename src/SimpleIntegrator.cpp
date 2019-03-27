#include <nori/integrator.h>
#include <nori/scene.h>
#include <minmax.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator
{
public:
	SimpleIntegrator(const PropertyList &props)
	{
		pointLightPosition = props.getPoint("position");
	 	pointLightEnergy = props.getColor("energy");
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override
	{
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		Point3f& hitPoint = its.p;
		Vector3f pointToLight = (pointLightPosition - hitPoint).normalized();
		
		Ray3f shadowRay = Ray3f(hitPoint, pointToLight);
		if (scene->rayIntersect(shadowRay))
			return Color3f(0.0f);

		float distanceSquare = pointToLight.dot(pointToLight);
		float cosTheta = max(0, its.shFrame.toLocal(pointToLight).normalized().z());

		return (pointLightEnergy / (4 * M_PI * M_PI)) * (cosTheta / distanceSquare);
	}

	std::string toString() const override
	{
		return "SimpleIntegrator";
	}

	Point3f pointLightPosition;
	Color3f pointLightEnergy;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END