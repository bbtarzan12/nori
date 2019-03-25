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

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <numeric>

NORI_NAMESPACE_BEGIN

constexpr uint32_t MAX_DEPTH = 10;
constexpr uint32_t MAX_CHILD = 8;

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build()
{
	if(m_octree)
		throw NoriException("Accel: only a single octree is supported!");

	m_octree = new Octree(m_mesh);
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const
{
	bool foundIntersection = false;  // Was an intersection found so far?
	uint32_t f = (uint32_t)-1;      // Triangle index of the closest intersection

	Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

	/* Brute force search through all triangles */
	//for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx)
	//{
	//	float u, v, t;
	//	if (m_mesh->rayIntersect(idx, ray, u, v, t))
	//	{
	//		/* An intersection was found! Can terminate
	//		   immediately if this is a shadow ray query */
	//		if (shadowRay)
	//			return true;
	//		ray.maxt = its.t = t;
	//		its.uv = Point2f(u, v);
	//		its.mesh = m_mesh;
	//		f = idx;
	//		foundIntersection = true;
	//	}
	//}

	uint32_t index;
	if (m_octree->search_tree(ray, index))
	{
		float u, v, t;
		if (m_mesh->rayIntersect(index, ray, u, v, t))
		{
			/* An intersection was found! Can terminate
			   immediately if this is a shadow ray query */
			if (shadowRay)
				return true;
			ray.maxt = its.t = t;
			its.uv = Point2f(u, v);
			its.mesh = m_mesh;
			f = index;
			foundIntersection = true;
		}
	}

	if (foundIntersection)
	{
		/* At this point, we now know that there is an intersection,
		   and we know the triangle index of the closest such intersection.

		   The following computes a number of additional properties which
		   characterize the intersection (normals, texture coordinates, etc..)
		*/

		/* Find the barycentric coordinates */
		Vector3f bary;
		bary << 1 - its.uv.sum(), its.uv;

		/* References to all relevant mesh buffers */
		const Mesh *mesh = its.mesh;
		const MatrixXf &V = mesh->getVertexPositions();
		const MatrixXf &N = mesh->getVertexNormals();
		const MatrixXf &UV = mesh->getVertexTexCoords();
		const MatrixXu &F = mesh->getIndices();

		/* Vertex indices of the triangle */
		uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

		Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

		/* Compute the intersection positon accurately
		   using barycentric coordinates */
		its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

		/* Compute proper texture coordinates if provided by the mesh */
		if (UV.size() > 0)
			its.uv = bary.x() * UV.col(idx0) +
			bary.y() * UV.col(idx1) +
			bary.z() * UV.col(idx2);

		/* Compute the geometry frame */
		its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());

		if (N.size() > 0)
		{
			/* Compute the shading frame. Note that for simplicity,
			   the current implementation doesn't attempt to provide
			   tangents that are continuous across the surface. That
			   means that this code will need to be modified to be able
			   use anisotropic BRDFs, which need tangent continuity */

			its.shFrame = Frame(
				(bary.x() * N.col(idx0) +
					bary.y() * N.col(idx1) +
					bary.z() * N.col(idx2)).normalized());
		}
		else
		{
			its.shFrame = its.geoFrame;
		}
	}

	return foundIntersection;
}

Octree::Octree(Mesh* mesh)
	: m_mesh(mesh), m_root(nullptr)
{
	if (!m_mesh)
		return;

	const BoundingBox3f& bbox = m_mesh->getBoundingBox();
	std::vector<uint32_t>* face_id = new std::vector<uint32_t>(m_mesh->getTriangleCount());
	std::iota(std::begin(*face_id), std::end(*face_id), 0);

	m_root = new Node(new BoundingBox3f(bbox.min, bbox.max), face_id, 0);
	m_root->generate_node(m_mesh);
}

bool nori::Octree::search_tree(const Ray3f& ray, uint32_t& index)
{
	float min_t = std::numeric_limits<float>::max();

	for (auto& child : *m_root->m_childs)
	{
		if (child != nullptr && child->m_bbox->rayIntersect(ray))
		{
			child->search_node(ray, m_mesh, index, min_t);
		}
	}

	return min_t != std::numeric_limits<float>::max();
}

void nori::Octree::Node::search_node(const Ray3f& ray, Mesh* mesh, uint32_t &index, float& min_t)
{
	if (is_leaf())
	{
		for (const auto& idx : *m_indicies)
		{
			float u, v, t;
			if (mesh->rayIntersect(idx, ray, u, v, t))
			{
				if (t < min_t)
				{
					index = idx;
					min_t = t;
				}
			}
		}
	}
	else
	{
		for (auto& child : *m_childs)
		{
			if (child != nullptr && child->m_bbox->rayIntersect(ray))
			{
				child->search_node(ray, mesh, index, min_t);
			}
		}
	}
}

void Octree::Node::generate_node(Mesh* mesh)
{
	if (m_depth >= MAX_DEPTH)
		return;

	m_childs = new std::vector<Node*>();

	const Point3f& center = m_bbox->getCenter();

	for (uint32_t idx = 0; idx < 8; idx++)
	{
		BoundingBox3f* child_bbox = new BoundingBox3f
		(
			Point3f
			(
				std::min(center[0], m_bbox->getCorner(idx)[0]),
				std::min(center[1], m_bbox->getCorner(idx)[1]),
				std::min(center[2], m_bbox->getCorner(idx)[2])
			),
			Point3f
			(
				std::max(center[0], m_bbox->getCorner(idx)[0]),
				std::max(center[1], m_bbox->getCorner(idx)[1]),
				std::max(center[2], m_bbox->getCorner(idx)[2])
			)
		);

		Node* child_node = new Node(child_bbox, new std::vector<uint32_t>(), m_depth + 1);
		m_childs->push_back(child_node);
	}

	for (const auto& index : *m_indicies)
	{
		for (auto& child : *m_childs)
		{
			if (child->m_bbox->overlaps(mesh->getBoundingBox(index)))
			{
				child->m_indicies->push_back(index);
			}
		}
	}

	delete m_indicies;
	m_indicies = nullptr;

	for (uint32_t idx = 0; idx < 8; idx++)
	{
		Node*& child_node = m_childs->at(idx);
		if (child_node->m_indicies->empty())
		{
			delete child_node;
			child_node = nullptr;
		}
		else if (child_node->m_indicies->size() > MAX_CHILD)
		{
			child_node->generate_node(mesh);
		}
	}
}

NORI_NAMESPACE_END
