/* mesh.c  -  mesh library  -  Public Domain  -  2020 Mattias Jansson
 *
 * This library provides a cross-platform mesh library in C11 providing
 * mesh handling functionality to applications base on our foundation library.
 *
 * The latest source code maintained by Mattias Jansson is always available at
 *
 * https://github.com/mjansson/mesh_lib
 *
 * This library is put in the public domain; you can redistribute it and/or modify it without any
 * restrictions.
 *
 */

#include "mesh.h"

#include <foundation/bucketarray.h>
#include <foundation/array.h>
#include <foundation/memory.h>
#include <foundation/math.h>
#include <foundation/log.h>
#include <foundation/stream.h>
#include <foundation/time.h>
#include <vector/vector.h>

static mesh_config_t _mesh_config;

int
mesh_module_initialize(mesh_config_t config) {
	_mesh_config = config;
	return 0;
}

void
mesh_module_finalize(void) {
}

bool
mesh_module_is_initialized(void) {
	return true;
}

void
mesh_module_parse_config(const char* path, size_t path_size, const char* buffer, size_t size,
                         const struct json_token_t* tokens, size_t tokens_count) {
	FOUNDATION_UNUSED(path);
	FOUNDATION_UNUSED(path_size);
	FOUNDATION_UNUSED(buffer);
	FOUNDATION_UNUSED(size);
	FOUNDATION_UNUSED(tokens);
	FOUNDATION_UNUSED(tokens_count);
}

mesh_t*
mesh_allocate(size_t expected_vertex_count, size_t expected_triangle_count) {
	mesh_t* mesh = memory_allocate(HASH_MESH, sizeof(mesh_t), 0, MEMORY_PERSISTENT);
	mesh_initialize(mesh, expected_vertex_count, expected_triangle_count);
	return mesh;
}

void
mesh_initialize(mesh_t* mesh, size_t expected_vertex_count, size_t expected_triangle_count) {
	bucketarray_initialize(&mesh->coordinate, sizeof(mesh_coordinate_t), expected_vertex_count / 4);
	bucketarray_initialize(&mesh->vertex, sizeof(mesh_vertex_t), expected_vertex_count / 4);
	bucketarray_initialize(&mesh->normal, sizeof(mesh_normal_t), expected_vertex_count / 8);
	bucketarray_initialize(&mesh->uv[0], sizeof(mesh_uv_t), expected_vertex_count / 8);
	bucketarray_initialize(&mesh->uv[1], sizeof(mesh_uv_t), expected_vertex_count / 8);
	bucketarray_initialize(&mesh->tangent, sizeof(mesh_tangent_t), expected_vertex_count / 8);
	bucketarray_initialize(&mesh->bitangent, sizeof(mesh_bitangent_t), expected_vertex_count / 8);
	bucketarray_initialize(&mesh->color, sizeof(mesh_color_t), expected_vertex_count / 8);
	bucketarray_initialize(&mesh->triangle, sizeof(mesh_triangle_t), expected_triangle_count / 4);
	mesh->attribute_vertex = nullptr;
	mesh->attribute_triangle = nullptr;
	mesh->partition = nullptr;
	mesh->topology = nullptr;
	mesh->next_lod = 0;
}

void
mesh_finalize(mesh_t* mesh) {
	mesh_t* first_mesh = mesh;
	while (mesh) {
		bucketarray_finalize(&mesh->coordinate);
		bucketarray_finalize(&mesh->vertex);
		bucketarray_finalize(&mesh->normal);
		bucketarray_finalize(&mesh->uv[0]);
		bucketarray_finalize(&mesh->uv[1]);
		bucketarray_finalize(&mesh->tangent);
		bucketarray_finalize(&mesh->bitangent);
		bucketarray_finalize(&mesh->color);
		bucketarray_finalize(&mesh->triangle);
		// TODO: free attributes
		// free mesh->attribute_vertex;
		// free mesh->attribute_triangle;
		if (mesh->partition) {
			bucketarray_finalize(&mesh->partition->coordinate_to_vertex);
			bucketarray_finalize(&mesh->partition->element);
		}
		memory_deallocate(mesh->partition);
		if (mesh->topology) {
			bucketarray_finalize(&mesh->topology->vertex_triangle_map);
			bucketarray_finalize(&mesh->topology->vertex_triangle_store);
		}
		memory_deallocate(mesh->topology);
		mesh_t* last_mesh = mesh;
		mesh = mesh->next_lod;
		if (last_mesh != first_mesh)
			memory_deallocate(last_mesh);
	}
}

void
mesh_deallocate(mesh_t* mesh) {
	if (mesh)
		mesh_finalize(mesh);
	memory_deallocate(mesh);
}

mesh_t*
mesh_clone(mesh_t* mesh) {
	mesh_t* clone = memory_allocate(HASH_MESH, sizeof(mesh_t), 0, MEMORY_PERSISTENT);
	bucketarray_initialize_copy(&clone->coordinate, &mesh->coordinate);
	bucketarray_initialize_copy(&clone->vertex, &mesh->vertex);
	bucketarray_initialize_copy(&clone->normal, &mesh->normal);
	bucketarray_initialize_copy(&clone->uv[0], &mesh->uv[0]);
	bucketarray_initialize_copy(&clone->uv[1], &mesh->uv[1]);
	bucketarray_initialize_copy(&clone->tangent, &mesh->tangent);
	bucketarray_initialize_copy(&clone->bitangent, &mesh->bitangent);
	bucketarray_initialize_copy(&clone->color, &mesh->color);
	bucketarray_initialize_copy(&clone->triangle, &mesh->triangle);
	clone->attribute_vertex = nullptr;
	clone->attribute_triangle = nullptr;
	clone->partition = nullptr;
	if (mesh->topology) {
		clone->topology = memory_allocate(HASH_MESH, sizeof(mesh_topology_t), 0, MEMORY_PERSISTENT);
		bucketarray_initialize_copy(&clone->topology->vertex_triangle_map, &mesh->topology->vertex_triangle_map);
		bucketarray_initialize_copy(&clone->topology->vertex_triangle_store, &mesh->topology->vertex_triangle_store);
	} else {
		clone->topology = nullptr;
	}
	clone->bounds_min = mesh->bounds_min;
	clone->bounds_max = mesh->bounds_max;
	clone->next_lod = 0;
	return clone;
}

void
mesh_set_coordinate_count(mesh_t* mesh, size_t count);

void
mesh_reserve_coordinate_count(mesh_t* mesh, size_t count);

void
mesh_set_normal_count(mesh_t* mesh, size_t count);

void
mesh_reserve_normal_count(mesh_t* mesh, size_t count);

void
mesh_set_vertex_count(mesh_t* mesh, size_t count);

void
mesh_reserve_vertex_count(mesh_t* mesh, size_t count);

void
mesh_set_triangle_count(mesh_t* mesh, size_t count);

void
mesh_reserve_triangle_count(mesh_t* mesh, size_t count);

void
mesh_calculate_bounds(mesh_t* mesh) {
	if (!mesh->coordinate.count)
		return;

	vector_t meshmin = vector(FLT_MAX, FLT_MAX, FLT_MAX, 1.0);
	vector_t meshmax = vector(-FLT_MAX, -FLT_MAX, -FLT_MAX, 1.0);

	for (size_t ivertex = 0; ivertex < mesh->vertex.count; ++ivertex) {
		mesh_vertex_t* vertex = bucketarray_get(&mesh->vertex, ivertex);
		mesh_coordinate_t* coordinate = bucketarray_get(&mesh->coordinate, vertex->coordinate);
		meshmin = vector_min(meshmin, *coordinate);
		meshmax = vector_max(meshmax, *coordinate);
	}

	mesh->bounds_min = meshmin;
	mesh->bounds_max = meshmax;
}

void
mesh_calculate_normals(mesh_t* mesh) {
	bool prev_normal_count = mesh->normal.count;

	// Clear any current normals
	bucketarray_resize(&mesh->normal, 0);
	if (prev_normal_count)
		bucketarray_resize_fill(&mesh->normal, prev_normal_count, 0);
	else
		bucketarray_resize_fill(&mesh->normal, mesh->coordinate.count, 0);

	vector_t minus_one = vector_uniform(-1.0f);
	vector_t one = vector_one();
	vector_t half = vector_half();

	// Generate triangle normals and add to vertex normals weighted by incident angle
	for (size_t itri = 0; itri < mesh->triangle.count; ++itri) {
		mesh_triangle_t* triangle = bucketarray_get(&mesh->triangle, itri);

		if (mesh_triangle_is_degenerate(mesh, triangle))
			continue;

		mesh_vertex_t* vertex[3] = {bucketarray_get(&mesh->vertex, triangle->vertex[0]),
		                            bucketarray_get(&mesh->vertex, triangle->vertex[1]),
		                            bucketarray_get(&mesh->vertex, triangle->vertex[2])};

		mesh_coordinate_t* coord[3] = {bucketarray_get(&mesh->coordinate, vertex[0]->coordinate),
		                               bucketarray_get(&mesh->coordinate, vertex[1]->coordinate),
		                               bucketarray_get(&mesh->coordinate, vertex[2]->coordinate)};

		// TODO: Support hard-angle cutoff
		if (!prev_normal_count) {
			vertex[0]->normal = vertex[0]->coordinate;
			vertex[1]->normal = vertex[1]->coordinate;
			vertex[2]->normal = vertex[2]->coordinate;
		}
		mesh_coordinate_t* normal[3] = {bucketarray_get(&mesh->normal, vertex[0]->normal),
		                                bucketarray_get(&mesh->normal, vertex[1]->normal),
		                                bucketarray_get(&mesh->normal, vertex[2]->normal)};

		vector_t edge0 = vector_sub(*coord[1], *coord[0]);
		vector_t edge1 = vector_sub(*coord[2], *coord[1]);
		vector_t edge2 = vector_sub(*coord[2], *coord[0]);

		vector_t triangle_normal = vector_normalize(vector_cross3(vector_normalize(edge0), vector_normalize(edge2)));

		edge0 = vector_normalize(edge0);
		edge1 = vector_normalize(edge1);
		edge2 = vector_normalize(edge2);

		vector_t e0_dot_e2 = vector_dot(edge0, edge2);
		vector_t e0_dot_e1 = vector_dot(vector_mul(edge0, minus_one), edge1);
		vector_t e1_dot_e2 = vector_dot(edge1, edge2);

		vector_t n0weighted = vector_mul(triangle_normal, vector_mul(vector_abs(vector_sub(one, e0_dot_e2)), half));
		vector_t n1weighted = vector_mul(triangle_normal, vector_mul(vector_abs(vector_sub(one, e0_dot_e1)), half));
		vector_t n2weighted = vector_mul(triangle_normal, vector_mul(vector_abs(vector_sub(one, e1_dot_e2)), half));

		*normal[0] = vector_add(*normal[0], n0weighted);
		*normal[1] = vector_add(*normal[1], n1weighted);
		*normal[2] = vector_add(*normal[2], n2weighted);
	}

	// Normalize
	for (size_t inormal = 0; inormal < mesh->normal.count; ++inormal) {
		mesh_normal_t* normal = bucketarray_get(&mesh->normal, inormal);
		vector_t sqrlen = vector_length_sqr(*normal);
		if (vector_x(sqrlen) > REAL_C(0.0001))
			*normal = vector_normalize(*normal);
		else
			*normal = vector(0, 1, 0, 0);
	}
}

void
mesh_calculate_tangents(mesh_t* mesh) {
	FOUNDATION_UNUSED(mesh);
}

static void
mesh_partition_insert(mesh_t* mesh, mesh_partition_t* partition, mesh_partition_element_t* element,
                      mesh_coordinate_t coordinate, unsigned int coordinate_index) {
	if (element->flags & MESH_PARTITION_FLAG_LEAF) {
		int iidx = 0;
		mesh_partition_leaf_t* leaf = &element->data.leaf;
		while ((iidx < MESH_PARTITION_NODE_MAX_ITEMS) && (leaf->coordinate[iidx] != MESH_INVALID_INDEX)) {
			if (leaf->coordinate[iidx] == coordinate_index)
				return;
			++iidx;
		}
		if (iidx < MESH_PARTITION_NODE_MAX_ITEMS) {
			// Leaf is not full
			leaf->coordinate[iidx] = coordinate_index;
		} else {
			// Leaf is full, split
			unsigned int child[8];
			vector_t split = vector_add(leaf->min, vector_scale(vector_sub(leaf->max, leaf->min), 0.5f));
			for (int ichild = 0; ichild < 8; ++ichild) {
				size_t element_count = partition->element.count;
				bucketarray_resize(&partition->element, element_count + 1);
				child[ichild] = (unsigned int)element_count;
				mesh_partition_element_t* child_element = bucketarray_get(&partition->element, element_count);
				child_element->flags = MESH_PARTITION_FLAG_LEAF;
				mesh_partition_leaf_t* child_leaf = &child_element->data.leaf;
				memset(child_leaf->coordinate, 0xFF, sizeof(child_leaf->coordinate));
				real min_x, min_y, min_z;
				real max_x, max_y, max_z;
				if (!(ichild & 1)) {
					min_x = vector_x(leaf->min);
					max_x = vector_x(split);
				} else {
					min_x = vector_x(split);
					max_x = vector_x(leaf->max);
				}
				if (!(ichild & 2)) {
					min_y = vector_y(leaf->min);
					max_y = vector_y(split);
				} else {
					min_y = vector_y(split);
					max_y = vector_y(leaf->max);
				}
				if (!(ichild & 4)) {
					min_z = vector_z(leaf->min);
					max_z = vector_z(split);
				} else {
					min_z = vector_z(split);
					max_z = vector_z(leaf->max);
				}
				child_leaf->min = vector(min_x, min_y, min_z, REAL_C(1.0));
				child_leaf->max = vector(max_x, max_y, max_z, REAL_C(1.0));
			}

			unsigned int leaf_coordinate[MESH_PARTITION_NODE_MAX_ITEMS];
			memcpy(leaf_coordinate, leaf->coordinate, sizeof(leaf_coordinate));

			element->flags &= ~(uint)MESH_PARTITION_FLAG_LEAF;
			mesh_partition_node_t* node = &element->data.node;
			node->split = split;
			memcpy(node->child, child, sizeof(node->child));

			for (iidx = 0; iidx < MESH_PARTITION_NODE_MAX_ITEMS; ++iidx) {
				unsigned int prev_coordinate_index = leaf_coordinate[iidx];
				mesh_coordinate_t* prev_coordinate = bucketarray_get(&mesh->coordinate, prev_coordinate_index);
				mesh_partition_insert(mesh, partition, element, *prev_coordinate, prev_coordinate_index);
			}
			mesh_partition_insert(mesh, partition, element, coordinate, coordinate_index);
		}
	} else {
		mesh_partition_node_t* node = &element->data.node;
		vectori_t less_than = vector_greater(coordinate, node->split);
		int32_t ichild = (vectori_x(less_than) & 1) | (vectori_y(less_than) & 2) | (vectori_z(less_than) & 4);
		mesh_partition_element_t* child_element = bucketarray_get(&partition->element, node->child[ichild]);
		mesh_partition_insert(mesh, partition, child_element, coordinate, coordinate_index);
	}
}

static void
mesh_validate_vertex_adjancency(mesh_t* mesh) {
	// Verify that all vertex rings are valid
	unsigned int* ring = nullptr;
	bool fail = false;
	for (unsigned int ivertex = 0; ivertex < mesh->vertex.count; ++ivertex) {
		mesh_vertex_t* vertex = bucketarray_get(&mesh->vertex, ivertex);
		array_clear(ring);
		while (!fail && (vertex->adjacent != ivertex)) {
			for (unsigned int iring = 0, ringcount = array_size(ring); iring < ringcount; ++iring) {
				if (ring[iring] == vertex->adjacent) {
					FOUNDATION_ASSERT_FAIL("Invalid vertex ring");
					fail = true;
					break;
				}
			}
			array_push(ring, vertex->adjacent);
			vertex = bucketarray_get(&mesh->vertex, vertex->adjacent);
		}
	}
	array_deallocate(ring);
}

static void
mesh_partition(mesh_t* mesh) {
	if (!mesh || !mesh->vertex.count || !mesh->coordinate.count || mesh->partition)
		return;

	vector_t bounds_dim = vector_sub(mesh->bounds_max, mesh->bounds_min);
	if ((vector_x(bounds_dim) <= 0) || (vector_y(bounds_dim) <= 0) || (vector_z(bounds_dim) <= 0)) {
		mesh_calculate_bounds(mesh);
		bounds_dim = vector_sub(mesh->bounds_max, mesh->bounds_min);
	}

	tick_t start = time_current();

	// Increase size by 5%
	vector_t offset = vector_scale(bounds_dim, 0.05f);
	vector_t bounds_min = vector_sub(mesh->bounds_min, offset);
	vector_t bounds_max = vector_add(mesh->bounds_max, offset);

	mesh_partition_t* partition = memory_allocate(HASH_MESH, sizeof(mesh_partition_t), 0, MEMORY_PERSISTENT);
	mesh->partition = partition;
	bucketarray_initialize(&partition->element, sizeof(mesh_partition_element_t),
	                       mesh->coordinate.count / (MESH_PARTITION_NODE_MAX_ITEMS * 2));

	partition->root.flags = MESH_PARTITION_FLAG_LEAF;
	partition->root.data.leaf.min = bounds_min;
	partition->root.data.leaf.max = bounds_max;
	memset(partition->root.data.leaf.coordinate, 0xFF, sizeof(partition->root.data.leaf.coordinate));
	// coordinate array now has MESH_INVALID_INDEX in all fields

	bucketarray_initialize(&mesh->partition->coordinate_to_vertex, sizeof(unsigned int),
	                       (size_t)1 << mesh->coordinate.bucket_shift);
	bucketarray_resize_fill(&mesh->partition->coordinate_to_vertex, mesh->coordinate.count, 0xFF);
	// coordinate to vertex array now has MESH_INVALID_INDEX in all fields

	for (unsigned int ivertex = 0; ivertex < mesh->vertex.count; ++ivertex) {
		mesh_vertex_t* vertex = bucketarray_get(&mesh->vertex, ivertex);
		mesh_coordinate_t* coordinate = bucketarray_get(&mesh->coordinate, vertex->coordinate);

		unsigned int* coordinate_vertex = bucketarray_get(&mesh->partition->coordinate_to_vertex, vertex->coordinate);
		if (*coordinate_vertex == MESH_INVALID_INDEX) {
			// First time this coordinate is used, insert into partitioning and link to vertex
			mesh_partition_insert(mesh, partition, &partition->root, *coordinate, vertex->coordinate);
			*coordinate_vertex = ivertex;
			vertex->adjacent = ivertex;
		} else {
			// Coordinate has been used before, add vertex to adjacent ring
			mesh_vertex_t* adjacent_vertex = bucketarray_get(&mesh->vertex, *coordinate_vertex);
			vertex->adjacent = adjacent_vertex->adjacent;
			adjacent_vertex->adjacent = ivertex;
		}
	}

	deltatime_t elapsed = time_elapsed(start);
	log_infof(HASH_MESH, STRING_CONST("Mesh partitioning created in %.2f seconds"), (double)elapsed);

	mesh_validate_vertex_adjancency(mesh);
}

static void
mesh_validate_topology(mesh_t* mesh) {
#if BUILD_ENABLE_ASSERT
	// Verify the topology
	for (unsigned int itri = 0; itri < mesh->triangle.count; ++itri) {
		mesh_triangle_t* triangle = bucketarray_get(&mesh->triangle, itri);
		for (unsigned int icorner = 0; icorner < 3; ++icorner) {
			unsigned int vertex_index = triangle->vertex[icorner];

			int found = 0;
			mesh_vertex_triangle_t* vertex_triangle =
			    bucketarray_get(&mesh->topology->vertex_triangle_map, vertex_index);
			size_t offset = vertex_triangle->offset;
			for (unsigned int it = 0, tcount = vertex_triangle->valence; it < tcount; ++it) {
				unsigned int* triangle_index = bucketarray_get(&mesh->topology->vertex_triangle_store, offset + it);
				if (*triangle_index == itri)
					++found;
			}
			FOUNDATION_ASSERT_MSG(found > 0,
			                      "Invalid vertex triangle map, triangle corner vertex does not map back to triangle");
			FOUNDATION_ASSERT_MSG(
			    found == 1, "Invalid vertex triangle map, triangle corner vertex map back to triangle multiple times");
		}
	}
	for (unsigned int ivertex = 0; ivertex < mesh->vertex.count; ++ivertex) {
		mesh_vertex_triangle_t* vertex_triangle = bucketarray_get(&mesh->topology->vertex_triangle_map, ivertex);
		size_t offset = vertex_triangle->offset;
		for (unsigned int it = 0, tcount = vertex_triangle->valence; it < tcount; ++it) {
			unsigned int* triangle_index = bucketarray_get(&mesh->topology->vertex_triangle_store, offset + it);
			mesh_triangle_t* triangle = bucketarray_get(&mesh->triangle, *triangle_index);
			bool found = false;
			for (unsigned int icorner = 0; icorner < 3; ++icorner) {
				unsigned int vertex_index = triangle->vertex[icorner];
				if (vertex_index == ivertex) {
					found = true;
					break;
				}
			}
			FOUNDATION_ASSERT_MSG(
			    found, "Invalid vertex triangle map, triangle list contains triangle that does not have vertex");
			for (unsigned int itnext = it + 1; itnext < tcount; ++itnext) {
				unsigned int* triangle_next_index =
				    bucketarray_get(&mesh->topology->vertex_triangle_store, offset + itnext);
				FOUNDATION_ASSERT_MSG(*triangle_index != *triangle_next_index,
				                      "Invalid vertex triangle map, duplicate entries");
			}
		}
	}
	for (unsigned int itri = 0; itri < mesh->triangle.count; ++itri) {
		mesh_triangle_t* triangle = bucketarray_get(&mesh->triangle, itri);
		for (unsigned int iedge = 0; iedge < 3; ++iedge) {
			if (triangle->adjacent[iedge] == MESH_INVALID_INDEX)
				continue;
			FOUNDATION_ASSERT_MSG(triangle->adjacent[iedge] < mesh->triangle.count,
			                      "Invalid adjacent triangle index, out of range");
			mesh_triangle_t* triangle_other = bucketarray_get(&mesh->triangle, triangle->adjacent[iedge]);
			bool found = false;
			for (unsigned int iotheredge = 0; iotheredge < 3; ++iotheredge) {
				if (triangle_other->adjacent[iotheredge] == itri)
					found = true;
			}
			FOUNDATION_ASSERT_MSG(
			    found, "Invalid triangle adjacency, adjacent triangle does not have origin triangle in adjacency list");
		}
	}
#else
	FOUNDATION_UNUSED(mesh);
#endif
}

static void
mesh_topology(mesh_t* mesh) {
	if (!mesh || !mesh->vertex.count || !mesh->coordinate.count)
		return;

	// Partition the mesh to get coordinate to vertex map and vertex adjacent rings setup
	if (!mesh->partition)
		mesh_partition(mesh);

	tick_t start = time_current();

	mesh_topology_t* topology =
	    memory_allocate(HASH_MESH, sizeof(mesh_topology_t), 0, MEMORY_PERSISTENT | MEMORY_ZERO_INITIALIZED);
	mesh->topology = topology;

	// First build the mapping from vertex to triangle
	bucketarray_initialize(&topology->vertex_triangle_map, sizeof(mesh_vertex_triangle_t), mesh->vertex.count / 4);
	bucketarray_initialize(&topology->vertex_triangle_store, sizeof(unsigned int), (mesh->triangle.count * 3) / 4);

	bucketarray_resize_fill(&topology->vertex_triangle_map, mesh->vertex.count, 0);

	// Collect the valence for all vertex, and reset triangle adjacency
	for (size_t itri = 0; itri < mesh->triangle.count; ++itri) {
		mesh_triangle_t* triangle = bucketarray_get(&mesh->triangle, itri);

		mesh_vertex_triangle_t* vertex[3] = {bucketarray_get(&topology->vertex_triangle_map, triangle->vertex[0]),
		                                     bucketarray_get(&topology->vertex_triangle_map, triangle->vertex[1]),
		                                     bucketarray_get(&topology->vertex_triangle_map, triangle->vertex[2])};

		++vertex[0]->valence;
		++vertex[1]->valence;
		++vertex[2]->valence;

		triangle->adjacent[0] = MESH_INVALID_INDEX;
		triangle->adjacent[1] = MESH_INVALID_INDEX;
		triangle->adjacent[2] = MESH_INVALID_INDEX;
	}

	// Allocate storage for a list of triangles for each vertex based on the valence
	size_t offset = 0;
	for (size_t ivertex = 0; ivertex < topology->vertex_triangle_map.count; ++ivertex) {
		mesh_vertex_triangle_t* vertex = bucketarray_get(&topology->vertex_triangle_map, ivertex);
		if (vertex->valence) {
			vertex->offset = offset;
			vertex->capacity = vertex->valence + 2;
			vertex->valence = 0;
			offset += vertex->capacity;
		}
	}

	bucketarray_resize(&topology->vertex_triangle_store, offset);

	// Collect the actual triangle list for each vertex
	for (unsigned int itri = 0; itri < mesh->triangle.count; ++itri) {
		mesh_triangle_t* triangle = bucketarray_get(&mesh->triangle, itri);

		mesh_vertex_triangle_t* vertex[3] = {bucketarray_get(&topology->vertex_triangle_map, triangle->vertex[0]),
		                                     bucketarray_get(&topology->vertex_triangle_map, triangle->vertex[1]),
		                                     bucketarray_get(&topology->vertex_triangle_map, triangle->vertex[2])};

		unsigned int* store = bucketarray_get(&topology->vertex_triangle_store, vertex[0]->offset + vertex[0]->valence);
		*store = itri;
		++vertex[0]->valence;

		store = bucketarray_get(&topology->vertex_triangle_store, vertex[1]->offset + vertex[1]->valence);
		*store = itri;
		++vertex[1]->valence;

		store = bucketarray_get(&topology->vertex_triangle_store, vertex[2]->offset + vertex[2]->valence);
		*store = itri;
		++vertex[2]->valence;
	}

	// Using the vertex to triangle mapping, collect adjacent triangles along each edge
	for (unsigned int itri = 0; itri < mesh->triangle.count; ++itri) {
		mesh_triangle_t* triangle = bucketarray_get(&mesh->triangle, itri);

		for (unsigned int iedge = 0; iedge < 3; ++iedge) {
			if (triangle->adjacent[iedge] != MESH_INVALID_INDEX)
				continue;

			// Triangle adjacency uses the coordinates, so get the coordinate index
			// for the start of the edge, and loop all vertex that share this coordinate
			unsigned int ifirstvertex = triangle->vertex[iedge];
			unsigned int istartvertex = ifirstvertex;
			mesh_vertex_t* vertex = bucketarray_get(&mesh->vertex, istartvertex);
			unsigned int istartcoordinate = vertex->coordinate;

			unsigned int inextvertex = triangle->vertex[(iedge + 1) % 3];
			mesh_vertex_t* next_vertex = bucketarray_get(&mesh->vertex, inextvertex);
			unsigned int inextcoordinate = next_vertex->coordinate;

			do {
				// Loop the triangles containing the current vertex and check if they have
				// an edge sharing the destination coordinate but with reversed winding
				mesh_vertex_triangle_t* vertex_triangle = bucketarray_get(&topology->vertex_triangle_map, istartvertex);
				vertex = bucketarray_get(&mesh->vertex, istartvertex);
				FOUNDATION_ASSERT(vertex->coordinate == istartcoordinate);

				unsigned int* other_triangle_index =
				    bucketarray_get(&topology->vertex_triangle_store, vertex_triangle->offset);
				for (unsigned int iotheridx = 0;
				     (iotheridx < vertex_triangle->valence) && (triangle->adjacent[iedge] == MESH_INVALID_INDEX);
				     ++iotheridx) {
					unsigned int iothertri = other_triangle_index[iotheridx];
					if (iothertri == itri)
						continue;

					// Check each edge in the other triangle if the vertex is the start vertex, and
					// the previous vertex has the same coordinate. If so, the edge runs between
					// the same coordinates as the edge currently processing, but with reversed winding
					mesh_triangle_t* other_triangle = bucketarray_get(&mesh->triangle, iothertri);
					for (unsigned int iotheredge = 0; iotheredge < 3; ++iotheredge) {
						unsigned int iothernextvertex = other_triangle->vertex[(iotheredge + 1) % 3];
						if (iothernextvertex != istartvertex)
							continue;

						unsigned int iotherstartvertex = other_triangle->vertex[iotheredge];
						mesh_vertex_t* other_start_vertex = bucketarray_get(&mesh->vertex, iotherstartvertex);
						if (other_start_vertex->coordinate != inextcoordinate)
							continue;

						if (other_triangle->adjacent[iotheredge] != MESH_INVALID_INDEX) {
							unsigned int ithirdtri = other_triangle->adjacent[iotheredge];
							mesh_triangle_t* third_triangle = bucketarray_get(&mesh->triangle, ithirdtri);

							// Pick the two triangles that have most aligned normals to be adjacent
							vector_t first_normal = mesh_triangle_normal(mesh, triangle);
							vector_t second_normal = mesh_triangle_normal(mesh, other_triangle);
							vector_t third_normal = mesh_triangle_normal(mesh, third_triangle);
							if (vectori_x(vector_greater(vector_dot(first_normal, second_normal),
							                             vector_dot(second_normal, third_normal)))) {
								// Replace previous adjacency
								triangle->adjacent[iedge] = iothertri;
								other_triangle->adjacent[iotheredge] = itri;

								// Reset adjacency of third triangle along the correct edge (winding same as original
								// triangle)
								unsigned int ithirdedge = 0;
								mesh_vertex_t* third_vertex[3] = {
								    bucketarray_get(&mesh->vertex, third_triangle->vertex[0]),
								    bucketarray_get(&mesh->vertex, third_triangle->vertex[1]),
								    bucketarray_get(&mesh->vertex, third_triangle->vertex[2])};
								if ((third_vertex[0]->coordinate == istartcoordinate) &&
								    (third_vertex[1]->coordinate == inextcoordinate)) {
									ithirdedge = 0;
								} else if ((third_vertex[1]->coordinate == istartcoordinate) &&
								           (third_vertex[2]->coordinate == inextcoordinate)) {
									ithirdedge = 1;
								} else if ((third_vertex[2]->coordinate == istartcoordinate) &&
								           (third_vertex[0]->coordinate == inextcoordinate)) {
									ithirdedge = 2;
								} else {
									FOUNDATION_ASSERT_FAIL("Internal triangle consistency failed");
									break;
								}
								third_triangle->adjacent[ithirdedge] = MESH_INVALID_INDEX;
								break;
							}
#if 0
									log_infof(
									    0,
									    STRING_CONST(
									        "Triangle %u vertex (%u, %u, %u) coord (%u, %u, %u) edge %u adjacent to:"),
									    itri, triangle->vertex[0], triangle->vertex[1], triangle->vertex[2],
									    bucketarray_get_as(mesh_vertex_t, &mesh->vertex, triangle->vertex[0])
									        ->coordinate,
									    bucketarray_get_as(mesh_vertex_t, &mesh->vertex, triangle->vertex[1])
									        ->coordinate,
									    bucketarray_get_as(mesh_vertex_t, &mesh->vertex, triangle->vertex[2])
									        ->coordinate,
									    iedge);
									log_infof(
									    0,
									    STRING_CONST("Triangle %u vertex (%u, %u, %u) coord (%u, %u, %u) edge %u "
									                 "already adjacent to:"),
									    iothertri, other_triangle->vertex[0], other_triangle->vertex[1],
									    other_triangle->vertex[2],
									    bucketarray_get_as(mesh_vertex_t, &mesh->vertex, other_triangle->vertex[0])
									        ->coordinate,
									    bucketarray_get_as(mesh_vertex_t, &mesh->vertex, other_triangle->vertex[1])
									        ->coordinate,
									    bucketarray_get_as(mesh_vertex_t, &mesh->vertex, other_triangle->vertex[2])
									        ->coordinate,
									    iotheredge);
									log_infof(
									    0, STRING_CONST("Triangle %u vertex (%u, %u, %u) coord (%u, %u, %u)"),
									    ithirdtri, third_triangle->vertex[0], third_triangle->vertex[1],
									    third_triangle->vertex[2],
									    bucketarray_get_as(mesh_vertex_t, &mesh->vertex, third_triangle->vertex[0])
									        ->coordinate,
									    bucketarray_get_as(mesh_vertex_t, &mesh->vertex, third_triangle->vertex[1])
									        ->coordinate,
									    bucketarray_get_as(mesh_vertex_t, &mesh->vertex, third_triangle->vertex[2])
									        ->coordinate);
									FOUNDATION_ASSERT_FAIL("Tringle edge split detected");
#endif
						} else {
							// We found an adjacent triangle along the edge, same coordinates
							// but reversed winding. Store adjacency for both triangles
							triangle->adjacent[iedge] = iothertri;
							other_triangle->adjacent[iotheredge] = itri;
							break;
						}
					}
				}

				istartvertex = vertex->adjacent;
			} while ((istartvertex != ifirstvertex) && (triangle->adjacent[iedge] == MESH_INVALID_INDEX));
		}
	}

	deltatime_t elapsed = time_elapsed(start);
	log_infof(HASH_MESH, STRING_CONST("Mesh topology created in %.2f seconds"), (double)elapsed);

	mesh_validate_topology(mesh);
}

void
mesh_merge_coordinate(mesh_t* mesh, real merge_coordinate_distance) {
	if (!mesh || !mesh->vertex.count || !mesh->coordinate.count)
		return;

	if (!mesh->partition)
		mesh_partition(mesh);
	FOUNDATION_UNUSED(merge_coordinate_distance);

	tick_t start = time_current();

	for (size_t inode = 0; inode < mesh->partition->element.count; ++inode) {
		mesh_partition_element_t* element = bucketarray_get(&mesh->partition->element, inode);
		if (!(element->flags & MESH_PARTITION_FLAG_LEAF))
			continue;

		mesh_partition_leaf_t* leaf = &element->data.leaf;
		for (unsigned int ifirstidx = 0; ifirstidx < MESH_PARTITION_NODE_MAX_ITEMS; ++ifirstidx) {
			unsigned int ifirst = leaf->coordinate[ifirstidx];
			if (ifirst == MESH_INVALID_INDEX)
				break;

			mesh_coordinate_t* first_coordinate = bucketarray_get(&mesh->coordinate, ifirst);

			for (unsigned int isecondidx = ifirstidx + 1; isecondidx < MESH_PARTITION_NODE_MAX_ITEMS;) {
				unsigned int isecond = leaf->coordinate[isecondidx];
				if (isecond == MESH_INVALID_INDEX)
					break;

				mesh_coordinate_t* second_coordinate = bucketarray_get(&mesh->coordinate, isecond);

				vector_t diff_sqr = vector_length_sqr(vector_sub(*first_coordinate, *second_coordinate));
				// if (math_real_eq(vector_x(diff_sqr), 0, 1000)) {
				if (vector_x(diff_sqr) < FLT_EPSILON) {
					// Coordinates are to be merged, update all vertex referring to the merged coordinate
					// and reset coordinate-to-vertex mapping for the merged coordinate
					unsigned int* start_vertex = bucketarray_get(&mesh->partition->coordinate_to_vertex, isecond);
					unsigned int ivertex = *start_vertex;
					do {
						mesh_vertex_t* vertex = bucketarray_get(&mesh->vertex, ivertex);
						vertex->coordinate = ifirst;
						ivertex = vertex->adjacent;
					} while (ivertex != *start_vertex);
					*start_vertex = MESH_INVALID_INDEX;

					// Splice in the vertex adjacent ring to the target coordinate vertex ring
					start_vertex = bucketarray_get(&mesh->partition->coordinate_to_vertex, ifirst);
					mesh_vertex_t* base_vertex = bucketarray_get(&mesh->vertex, *start_vertex);

					mesh_vertex_t* splice_vertex = bucketarray_get(&mesh->vertex, ivertex);
					unsigned int inextvertex = base_vertex->adjacent;
					base_vertex->adjacent = splice_vertex->adjacent;
					splice_vertex->adjacent = inextvertex;

					// Remove the coordinate from the partitioning by swap-with-last
					unsigned int ilastidx = MESH_PARTITION_NODE_MAX_ITEMS - 1;
					for (; ilastidx > isecondidx; --ilastidx) {
						if (leaf->coordinate[ilastidx] != MESH_INVALID_INDEX)
							break;
					}
					leaf->coordinate[isecondidx] = leaf->coordinate[ilastidx];
					leaf->coordinate[ilastidx] = MESH_INVALID_INDEX;
				} else {
					++isecondidx;
				}
			}
		}
	}

	deltatime_t elapsed = time_elapsed(start);
	log_infof(HASH_MESH, STRING_CONST("Mesh coordinates merged in %.2f seconds"), (double)elapsed);
}

void
mesh_merge_vertex(mesh_t* mesh) {
	if (!mesh)
		return;

	if (!mesh->partition)
		mesh_partition(mesh);

	if (!mesh->topology)
		mesh_topology(mesh);

	tick_t start = time_current();

	vector_t tolerance = vector_uniform(FLT_EPSILON);
	for (unsigned int ivertex = 0; ivertex < mesh->vertex.count; ++ivertex) {
		mesh_vertex_t* vertex = bucketarray_get(&mesh->vertex, ivertex);
		unsigned int iadjacent = vertex->adjacent;
		while (iadjacent != ivertex) {
			mesh_vertex_t* vertex_adjacent = bucketarray_get(&mesh->vertex, iadjacent);

			int32_t equal = 1;
			if (mesh->normal.count && (vertex->normal != vertex_adjacent->normal)) {
				mesh_normal_t* normal_first = bucketarray_get(&mesh->normal, vertex->normal);
				mesh_normal_t* normal_second = bucketarray_get(&mesh->normal, vertex_adjacent->normal);
				vector_t diff = vector_sub(*normal_first, *normal_second);
				if (vectori_x(vector_greater(vector_length_sqr(diff), tolerance)))
					equal = 0;
			}
			if (equal && mesh->uv[0].count && (vertex->uv[0] != vertex_adjacent->uv[0])) {
				mesh_uv_t* uv_first = bucketarray_get(&mesh->uv[0], vertex->uv[0]);
				mesh_uv_t* uv_second = bucketarray_get(&mesh->uv[0], vertex_adjacent->uv[0]);
				if (!math_real_eq(uv_first->u, uv_second->u, 1000) || !math_real_eq(uv_first->v, uv_second->v, 1000))
					equal = 0;
			}
			if (equal && mesh->uv[1].count && (vertex->uv[1] != vertex_adjacent->uv[1])) {
				mesh_uv_t* uv_first = bucketarray_get(&mesh->uv[1], vertex->uv[1]);
				mesh_uv_t* uv_second = bucketarray_get(&mesh->uv[1], vertex_adjacent->uv[1]);
				if (!math_real_eq(uv_first->u, uv_second->u, 1000) || !math_real_eq(uv_first->v, uv_second->v, 1000))
					equal = 0;
			}

			if (equal) {
				// Merge adjacent vertex into this vertex, update triangles and topology
				mesh_vertex_triangle_t* adjacent_vertex_triangle =
				    bucketarray_get(&mesh->topology->vertex_triangle_map, iadjacent);
				unsigned int* adjacent_triangle_index =
				    bucketarray_get(&mesh->topology->vertex_triangle_store, adjacent_vertex_triangle->offset);
				for (unsigned int itri = 0; itri < adjacent_vertex_triangle->valence; ++itri) {
					mesh_triangle_t* triangle = bucketarray_get(&mesh->triangle, adjacent_triangle_index[itri]);
					if (triangle->vertex[0] == iadjacent)
						triangle->vertex[0] = ivertex;
					if (triangle->vertex[1] == iadjacent)
						triangle->vertex[1] = ivertex;
					if (triangle->vertex[2] == iadjacent)
						triangle->vertex[2] = ivertex;
				}

				mesh_vertex_triangle_t* vertex_triangle =
				    bucketarray_get(&mesh->topology->vertex_triangle_map, ivertex);
				unsigned int total_valence = vertex_triangle->valence + adjacent_vertex_triangle->valence;
				if (total_valence > vertex_triangle->capacity) {
					size_t new_offset = mesh->topology->vertex_triangle_store.count;
					total_valence += 2;
					bucketarray_resize(&mesh->topology->vertex_triangle_store, new_offset + total_valence);
					memcpy(bucketarray_get(&mesh->topology->vertex_triangle_store, new_offset),
					       bucketarray_get(&mesh->topology->vertex_triangle_store, vertex_triangle->offset),
					       sizeof(unsigned int) * vertex_triangle->valence);
					vertex_triangle->offset = new_offset;
					vertex_triangle->capacity = total_valence;
				}
				unsigned int* triangle_index =
				    bucketarray_get(&mesh->topology->vertex_triangle_store, vertex_triangle->offset);
				memcpy(triangle_index + vertex_triangle->valence, adjacent_triangle_index,
				       sizeof(unsigned int) * adjacent_vertex_triangle->valence);
				vertex_triangle->valence += adjacent_vertex_triangle->valence;
				adjacent_vertex_triangle->valence = 0;
			}

			iadjacent = vertex_adjacent->adjacent;
		}
	}

	deltatime_t elapsed = time_elapsed(start);
	log_infof(HASH_MESH, STRING_CONST("Mesh vertex merged in %.2f seconds"), (double)elapsed);
}

void
mesh_compact(mesh_t* mesh) {
	if (!mesh)
		return;

	if (!mesh->triangle.count || !mesh->vertex.count)
		return;

	// Partition the mesh to get coordinate to vertex map and vertex adjacent rings setup
	if (!mesh->partition)
		mesh_partition(mesh);

	// Calculate mesh topology to get vertex to triangle mappings
	if (!mesh->topology)
		mesh_topology(mesh);

	tick_t start = time_current();

	// Discard unused data, will be recreated
	bucketarray_finalize(&mesh->topology->vertex_triangle_store);

	// First compact triangle storage and create a map from old triangle index
	// to new triangle index, cleaning out degenerate triangles
	bucketarray_t triangle_map;
	bucketarray_initialize(&triangle_map, sizeof(unsigned int), mesh->triangle.count / 4);
	bucketarray_resize(&triangle_map, mesh->triangle.count);

	bucketarray_t old_triangle_array = mesh->triangle;
	bucketarray_initialize(&mesh->triangle, sizeof(mesh_triangle_t), old_triangle_array.count / 4);

	for (unsigned int ioldtri = 0; ioldtri < old_triangle_array.count; ++ioldtri) {
		unsigned int* mapping = bucketarray_get(&triangle_map, ioldtri);
		mesh_triangle_t* triangle = bucketarray_get(&old_triangle_array, ioldtri);
		bool valid = !mesh_triangle_is_degenerate(mesh, triangle);
		if (valid) {
			*mapping = (unsigned int)mesh->triangle.count;
			bucketarray_push(&mesh->triangle, triangle);
		} else {
			// When removing a triangle, reduce valence of affected vertex. Do not bother
			// with updating triangle list yet, will be done when compacting storage
			*mapping = MESH_INVALID_INDEX;
			for (unsigned int icorner = 0; icorner < 3; ++icorner) {
				if (triangle->vertex[icorner] == MESH_INVALID_INDEX)
					continue;
				mesh_vertex_triangle_t* vertex_triangle =
				    bucketarray_get(&mesh->topology->vertex_triangle_map, triangle->vertex[icorner]);
				FOUNDATION_ASSERT(vertex_triangle->valence);
				--vertex_triangle->valence;
			}
		}
	}

	bucketarray_finalize(&old_triangle_array);
	log_infof(HASH_MESH, STRING_CONST("Mesh compacted triangle array: %" PRIsize " to %" PRIsize),
	          old_triangle_array.count, mesh->triangle.count);

	// Now compact vertex storage and create a map from old vertex index to new
	// vertex index, cleaning out unused vertex and calculating storage
	bucketarray_t vertex_map;
	bucketarray_initialize(&vertex_map, sizeof(unsigned int), mesh->vertex.count / 4);
	bucketarray_resize(&vertex_map, mesh->vertex.count);

	size_t new_vertex_count = 0;
	for (unsigned int ioldvertex = 0; ioldvertex < mesh->vertex.count; ++ioldvertex) {
		mesh_vertex_triangle_t* vertex_triangle = bucketarray_get(&mesh->topology->vertex_triangle_map, ioldvertex);
		if (vertex_triangle->valence > 0)
			++new_vertex_count;
	}

	bucketarray_t old_vertex_array = mesh->vertex;
	bucketarray_initialize(&mesh->vertex, sizeof(mesh_vertex_t), new_vertex_count / 4);

	bucketarray_t old_vertex_triangle_map = mesh->topology->vertex_triangle_map;
	bucketarray_initialize(&mesh->topology->vertex_triangle_map, sizeof(mesh_vertex_triangle_t), new_vertex_count / 4);

	size_t triangle_map_capacity = 0;
	for (unsigned int ioldvertex = 0; ioldvertex < old_vertex_array.count; ++ioldvertex) {
		unsigned int* mapping = bucketarray_get(&vertex_map, ioldvertex);
		mesh_vertex_t* vertex = bucketarray_get(&old_vertex_array, ioldvertex);
		mesh_vertex_triangle_t* vertex_triangle = bucketarray_get(&old_vertex_triangle_map, ioldvertex);
		if (vertex_triangle->valence > 0) {
			*mapping = (unsigned int)mesh->vertex.count;
			vertex_triangle->capacity = vertex_triangle->valence;
			vertex_triangle->offset = triangle_map_capacity;
			vertex_triangle->valence = 0;
			triangle_map_capacity += vertex_triangle->capacity;
			bucketarray_push(&mesh->vertex, vertex);
			bucketarray_push(&mesh->topology->vertex_triangle_map, vertex_triangle);
		} else {
			*mapping = MESH_INVALID_INDEX;
		}
	}

	bucketarray_finalize(&old_vertex_triangle_map);
	log_infof(HASH_MESH, STRING_CONST("Mesh compacted vertex array: %" PRIsize " to %" PRIsize), old_vertex_array.count,
	          mesh->vertex.count);

	// Update all vertex adjacent links, removing any discarded vertex
	for (unsigned int ivertex = 0; ivertex < mesh->vertex.count; ++ivertex) {
		mesh_vertex_t* vertex = bucketarray_get(&mesh->vertex, ivertex);
		unsigned int first_adjacent = vertex->adjacent;
		unsigned int adjacent = first_adjacent;
		unsigned int valid = 0;
		do {
			unsigned int* mapped_adjacent = bucketarray_get(&vertex_map, adjacent);
			if (*mapped_adjacent != MESH_INVALID_INDEX) {
				adjacent = *mapped_adjacent;
				valid = 1;
				break;
			}
			mesh_vertex_t* old_vertex = bucketarray_get(&old_vertex_array, adjacent);
			adjacent = old_vertex->adjacent;
		} while (adjacent != first_adjacent);
		if (!valid)
			adjacent = ivertex;
		vertex->adjacent = adjacent;
	}

	// Update all triangles with new vertex index and new adjacent triangle index
	// Update all vertex to triangle maps with new triangle index and compact storage,
	// recreate triangle list for each vertex
	bucketarray_initialize(&mesh->topology->vertex_triangle_store, sizeof(unsigned int), triangle_map_capacity / 4);
	bucketarray_resize(&mesh->topology->vertex_triangle_store, triangle_map_capacity);

	for (unsigned int itri = 0; itri < mesh->triangle.count; ++itri) {
		mesh_triangle_t* triangle = bucketarray_get(&mesh->triangle, itri);
		for (unsigned int icorner = 0; icorner < 3; ++icorner) {
			unsigned int* vertex_mapping = bucketarray_get(&vertex_map, triangle->vertex[icorner]);
			unsigned int ivertex = *vertex_mapping;
			FOUNDATION_ASSERT(ivertex != MESH_INVALID_INDEX);
			triangle->vertex[icorner] = ivertex;

			mesh_vertex_triangle_t* vertex_triangle = bucketarray_get(&mesh->topology->vertex_triangle_map, ivertex);
			unsigned int* vertex_triangle_index = bucketarray_get(&mesh->topology->vertex_triangle_store,
			                                                      vertex_triangle->offset + vertex_triangle->valence);
			FOUNDATION_ASSERT(vertex_triangle->valence < vertex_triangle->capacity);
			*vertex_triangle_index = itri;
			++vertex_triangle->valence;

			unsigned int iedge = icorner;
			if (triangle->adjacent[iedge] != MESH_INVALID_INDEX) {
				unsigned int* triangle_mapping = bucketarray_get(&triangle_map, triangle->adjacent[iedge]);
				triangle->adjacent[iedge] = *triangle_mapping;
			}
		}
	}

	bucketarray_finalize(&triangle_map);

	// Compact coordinate array and update coordinate to vertex map
	size_t new_coordinate_count = 0;
	for (unsigned int icoord = 0; icoord < mesh->partition->coordinate_to_vertex.count; ++icoord) {
		unsigned int* coordinate_vertex = bucketarray_get(&mesh->partition->coordinate_to_vertex, icoord);
		unsigned int ivertex = *coordinate_vertex;
		if (ivertex != MESH_INVALID_INDEX) {
			bool mapped = false;
			do {
				unsigned int* mapped_vertex = bucketarray_get(&vertex_map, ivertex);
				if (*mapped_vertex != MESH_INVALID_INDEX) {
					ivertex = *mapped_vertex;
					mapped = true;
					FOUNDATION_ASSERT(ivertex < mesh->vertex.count);
					break;
				}
				mesh_vertex_t* old_vertex = bucketarray_get(&old_vertex_array, ivertex);
				ivertex = old_vertex->adjacent;
			} while ((ivertex != *coordinate_vertex) && (ivertex != MESH_INVALID_INDEX));
			if (!mapped)
				ivertex = MESH_INVALID_INDEX;
		}
		FOUNDATION_ASSERT((ivertex < mesh->vertex.count) || (ivertex == MESH_INVALID_INDEX));
		*coordinate_vertex = ivertex;
		if (ivertex != MESH_INVALID_INDEX)
			++new_coordinate_count;
	}

	bucketarray_t old_coordinate_array = mesh->coordinate;
	bucketarray_initialize(&mesh->coordinate, sizeof(mesh_coordinate_t), new_coordinate_count / 4);

	bucketarray_t old_coordinate_to_vertex_array = mesh->partition->coordinate_to_vertex;
	bucketarray_initialize(&mesh->partition->coordinate_to_vertex, sizeof(unsigned int), new_coordinate_count / 4);

	for (unsigned int icoord = 0; icoord < old_coordinate_to_vertex_array.count; ++icoord) {
		unsigned int* coordinate_vertex = bucketarray_get(&old_coordinate_to_vertex_array, icoord);
		if (*coordinate_vertex != MESH_INVALID_INDEX) {
			unsigned int icoordinate = (unsigned int)mesh->coordinate.count;
			unsigned int ifirstvertex = *coordinate_vertex;
			unsigned int ivertex = ifirstvertex;

			bucketarray_push(&mesh->coordinate, bucketarray_get(&old_coordinate_array, icoord));
			bucketarray_push(&mesh->partition->coordinate_to_vertex, coordinate_vertex);

			do {
				mesh_vertex_t* vertex = bucketarray_get(&mesh->vertex, ivertex);
				vertex->coordinate = icoordinate;
				ivertex = vertex->adjacent;
			} while (ivertex != ifirstvertex);
		}
	}

	bucketarray_finalize(&old_coordinate_array);
	bucketarray_finalize(&old_coordinate_to_vertex_array);
	log_infof(HASH_MESH, STRING_CONST("Mesh compacted coordinate array: %" PRIsize " to %" PRIsize),
	          old_coordinate_array.count, mesh->coordinate.count);

	bucketarray_finalize(&old_vertex_array);
	bucketarray_finalize(&vertex_map);

	// TODO: attributes
	if (mesh->normal.count) {
	}

	for (unsigned int ilayer = 0; ilayer < 2; ++ilayer) {
		if (mesh->uv[ilayer].count) {
		}
	}

	deltatime_t elapsed = time_elapsed(start);
	log_infof(HASH_MESH, STRING_CONST("Mesh compacted in %.2f seconds"), (double)elapsed);

	mesh_validate_vertex_adjancency(mesh);
	mesh_validate_topology(mesh);
}

void
mesh_calculate_topology(mesh_t* mesh) {
	if (mesh->topology) {
		bucketarray_finalize(&mesh->topology->vertex_triangle_map);
		bucketarray_finalize(&mesh->topology->vertex_triangle_store);
		memory_deallocate(mesh->topology);
		mesh->topology = nullptr;
	}

	mesh_topology(mesh);
}

int
mesh_valid(mesh_t* mesh) {
	FOUNDATION_UNUSED(mesh);
	return 0;
}

static void
mesh_dump_bucketarray(stream_t* stream, bucketarray_t* array) {
	uint64_t count = array->count;
	stream_write_uint64(stream, count);

	size_t bucket_element_count = (size_t)1 << array->bucket_shift;
	for (size_t ibucket = 0; count > 0; ++ibucket) {
		if (count >= bucket_element_count) {
			stream_write(stream, array->bucket[ibucket], array->element_size * bucket_element_count);
			count -= bucket_element_count;
		} else {
			stream_write(stream, array->bucket[ibucket], array->element_size * count);
			count = 0;
		}
	}
}

static void
mesh_undump_bucketarray(stream_t* stream, bucketarray_t* array) {
	uint64_t count = stream_read_uint64(stream);
	bucketarray_initialize(array, array->element_size, (size_t)count / 4);
	bucketarray_resize(array, count);

	size_t bucket_element_count = (size_t)1 << array->bucket_shift;
	for (size_t ibucket = 0; count > 0; ++ibucket) {
		if (count >= bucket_element_count) {
			stream_read(stream, array->bucket[ibucket], array->element_size * bucket_element_count);
			count -= bucket_element_count;
		} else {
			stream_read(stream, array->bucket[ibucket], array->element_size * count);
			count = 0;
		}
	}
}

void
mesh_dump(mesh_t* mesh, stream_t* stream) {
	stream_set_binary(stream, true);

	mesh_dump_bucketarray(stream, &mesh->coordinate);
	mesh_dump_bucketarray(stream, &mesh->normal);
	mesh_dump_bucketarray(stream, &mesh->uv[0]);
	mesh_dump_bucketarray(stream, &mesh->uv[1]);
	mesh_dump_bucketarray(stream, &mesh->tangent);
	mesh_dump_bucketarray(stream, &mesh->bitangent);
	mesh_dump_bucketarray(stream, &mesh->vertex);
	mesh_dump_bucketarray(stream, &mesh->triangle);

	stream_write_bool(stream, mesh->partition ? true : false);
	if (mesh->partition) {
		stream_write(stream, &mesh->partition->root, sizeof(mesh_partition_element_t));
		mesh_dump_bucketarray(stream, &mesh->partition->element);
		mesh_dump_bucketarray(stream, &mesh->partition->coordinate_to_vertex);
	}

	stream_write_bool(stream, mesh->topology ? true : false);
	if (mesh->topology) {
		mesh_dump_bucketarray(stream, &mesh->topology->vertex_triangle_map);
		mesh_dump_bucketarray(stream, &mesh->topology->vertex_triangle_store);
	}

	stream_write(stream, &mesh->bounds_min, sizeof(mesh->bounds_min));
	stream_write(stream, &mesh->bounds_max, sizeof(mesh->bounds_max));
}

mesh_t*
mesh_undump(stream_t* stream) {
	stream_set_binary(stream, true);

	mesh_t* mesh = mesh_allocate(0, 0);
	mesh_undump_bucketarray(stream, &mesh->coordinate);
	mesh_undump_bucketarray(stream, &mesh->normal);
	mesh_undump_bucketarray(stream, &mesh->uv[0]);
	mesh_undump_bucketarray(stream, &mesh->uv[1]);
	mesh_undump_bucketarray(stream, &mesh->tangent);
	mesh_undump_bucketarray(stream, &mesh->bitangent);
	mesh_undump_bucketarray(stream, &mesh->vertex);
	mesh_undump_bucketarray(stream, &mesh->triangle);

	if (stream_read_bool(stream)) {
		mesh->partition = memory_allocate(HASH_MESH, sizeof(mesh_partition_t), 0, MEMORY_PERSISTENT);
		bucketarray_initialize(&mesh->partition->element, sizeof(mesh_partition_element_t), 0);
		bucketarray_initialize(&mesh->partition->coordinate_to_vertex, sizeof(unsigned int), 0);

		stream_read(stream, &mesh->partition->root, sizeof(mesh_partition_element_t));
		mesh_undump_bucketarray(stream, &mesh->partition->element);
		mesh_undump_bucketarray(stream, &mesh->partition->coordinate_to_vertex);
	}

	if (stream_read_bool(stream)) {
		mesh->topology =
		    memory_allocate(HASH_MESH, sizeof(mesh_topology_t), 0, MEMORY_PERSISTENT | MEMORY_ZERO_INITIALIZED);
		bucketarray_initialize(&mesh->topology->vertex_triangle_map, sizeof(mesh_vertex_triangle_t), 0);
		bucketarray_initialize(&mesh->topology->vertex_triangle_store, sizeof(unsigned int), 0);

		mesh_undump_bucketarray(stream, &mesh->topology->vertex_triangle_map);
		mesh_undump_bucketarray(stream, &mesh->topology->vertex_triangle_store);
	}

	stream_read(stream, &mesh->bounds_min, sizeof(mesh->bounds_min));
	stream_read(stream, &mesh->bounds_max, sizeof(mesh->bounds_max));

	return mesh;
}

bool
mesh_triangle_is_degenerate(mesh_t* mesh, mesh_triangle_t* triangle) {
	if ((triangle->vertex[0] == triangle->vertex[1]) || (triangle->vertex[0] == triangle->vertex[2]) ||
	    (triangle->vertex[1] == triangle->vertex[2]))
		return true;

	mesh_vertex_t* vertex[3] = {bucketarray_get(&mesh->vertex, triangle->vertex[0]),
	                            bucketarray_get(&mesh->vertex, triangle->vertex[1]),
	                            bucketarray_get(&mesh->vertex, triangle->vertex[2])};
	if ((vertex[0]->coordinate == vertex[1]->coordinate) || (vertex[0]->coordinate == vertex[2]->coordinate) ||
	    (vertex[1]->coordinate == vertex[2]->coordinate))
		return true;
	return false;
}

vector_t
mesh_triangle_normal(mesh_t* mesh, mesh_triangle_t* triangle) {
	mesh_vertex_t* vertex[3] = {bucketarray_get(&mesh->vertex, triangle->vertex[0]),
	                            bucketarray_get(&mesh->vertex, triangle->vertex[1]),
	                            bucketarray_get(&mesh->vertex, triangle->vertex[2])};
	if ((vertex[0]->coordinate == vertex[1]->coordinate) || (vertex[0]->coordinate == vertex[2]->coordinate) ||
	    (vertex[1]->coordinate == vertex[2]->coordinate))
		return vector_zero();

	mesh_coordinate_t* coordinate[3] = {bucketarray_get(&mesh->coordinate, vertex[0]->coordinate),
	                                    bucketarray_get(&mesh->coordinate, vertex[1]->coordinate),
	                                    bucketarray_get(&mesh->coordinate, vertex[2]->coordinate)};
	vector_t edge0 = vector_sub(*coordinate[1], *coordinate[0]);
	vector_t edge1 = vector_sub(*coordinate[2], *coordinate[0]);
	return vector_normalize(vector_cross3(edge0, edge1));
}

void
mesh_merge_mesh(mesh_t* mesh, mesh_t* additional) {
	if (!mesh || !additional)
		return;

	// TODO: Merge all attributes
	size_t coordinate_offset = mesh->coordinate.count;
	bucketarray_append(&mesh->coordinate, &additional->coordinate);

	bool additional_normals = (additional->normal.count > 0);
	size_t normal_offset = mesh->normal.count;
	if (additional_normals)
		bucketarray_append(&mesh->normal, &additional->normal);

	bool additional_tangents = (additional->tangent.count > 0);
	size_t tangent_offset = mesh->tangent.count;
	if (additional_tangents)
		bucketarray_append(&mesh->tangent, &additional->tangent);

	bool additional_bitangents = (additional->bitangent.count > 0);
	size_t bitangent_offset = mesh->bitangent.count;
	if (additional_bitangents)
		bucketarray_append(&mesh->bitangent, &additional->bitangent);

	bool additional_uvs[2] = {additional->uv[0].count > 0, additional->uv[1].count > 0};
	size_t uv_offset[2] = {mesh->uv[0].count, mesh->uv[1].count};
	if (additional_uvs[0])
		bucketarray_append(&mesh->uv[0], &additional->uv[0]);
	if (additional_uvs[1])
		bucketarray_append(&mesh->uv[1], &additional->uv[1]);

	size_t vertex_offset = mesh->vertex.count;
	bucketarray_append(&mesh->vertex, &additional->vertex);
	for (size_t ivert = vertex_offset; ivert < mesh->vertex.count; ++ivert) {
		mesh_vertex_t* vertex = bucketarray_get_as(mesh_vertex_t, &mesh->vertex, ivert);
		vertex->coordinate += (uint)coordinate_offset;
		if (additional_normals)
			vertex->normal += (uint)normal_offset;
		if (additional_tangents)
			vertex->tangent += (uint)tangent_offset;
		if (additional_bitangents)
			vertex->bitangent += (uint)bitangent_offset;
		if (additional_uvs[0])
			vertex->uv[0] += (uint)uv_offset[0];
		if (additional_uvs[0])
			vertex->uv[1] += (uint)uv_offset[1];
	}

	size_t triangle_offset = mesh->triangle.count;
	bucketarray_append(&mesh->triangle, &additional->triangle);
	for (size_t itri = triangle_offset; itri < mesh->triangle.count; ++itri) {
		mesh_triangle_t* triangle = bucketarray_get_as(mesh_triangle_t, &mesh->triangle, itri);
		triangle->vertex[0] += (uint)vertex_offset;
		triangle->vertex[1] += (uint)vertex_offset;
		triangle->vertex[2] += (uint)vertex_offset;
	}

	mesh->bounds_max = vector_max(mesh->bounds_max, additional->bounds_max);
	mesh->bounds_min = vector_min(mesh->bounds_min, additional->bounds_min);
}

void
mesh_merge_transformed_mesh(mesh_t* mesh, mesh_t* additional, matrix_t transform) {
	if (!mesh || !additional)
		return;

	vector_t additional_bounds_min = vector(REAL_MAX, REAL_MAX, REAL_MAX, REAL_MAX);
	vector_t additional_bounds_max = vector(-REAL_MAX, -REAL_MAX, -REAL_MAX, -REAL_MAX);

	// TODO: Merge attributes
	size_t coordinate_offset = mesh->coordinate.count;
	for (size_t icoord = 0; icoord < additional->coordinate.count; ++icoord) {
		vector_t coordinate = *bucketarray_get_as(vector_t, &additional->coordinate, icoord);
		coordinate = vector_transform(coordinate, transform);
		bucketarray_push(&mesh->coordinate, &coordinate);
	}

	size_t vertex_offset = mesh->vertex.count;
	bucketarray_append(&mesh->vertex, &additional->vertex);
	for (size_t ivert = vertex_offset; ivert < mesh->vertex.count; ++ivert)
		bucketarray_get_as(mesh_vertex_t, &mesh->vertex, ivert)->coordinate += (uint)coordinate_offset;

	size_t triangle_offset = mesh->triangle.count;
	bucketarray_append(&mesh->triangle, &additional->triangle);
	for (size_t itri = triangle_offset; itri < mesh->triangle.count; ++itri) {
		mesh_triangle_t* triangle = bucketarray_get_as(mesh_triangle_t, &mesh->triangle, itri);
		triangle->vertex[0] += (uint)vertex_offset;
		triangle->vertex[1] += (uint)vertex_offset;
		triangle->vertex[2] += (uint)vertex_offset;
	}

	mesh->bounds_max = vector_max(mesh->bounds_max, additional_bounds_max);
	mesh->bounds_min = vector_min(mesh->bounds_min, additional_bounds_min);
}
