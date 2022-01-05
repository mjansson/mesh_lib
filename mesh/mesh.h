/* mesh.h  -  mesh library  -  Public Domain  -  2020 Mattias Jansson
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

#pragma once

/*! \file mesh.h
    Mesh library entry points */

#include <mesh/types.h>
#include <mesh/hashstrings.h>

/*! Initialize mesh library
    \return 0 if success, <0 if error */
MESH_API int
mesh_module_initialize(mesh_config_t config);

/*! Finalize mesh library */
MESH_API void
mesh_module_finalize(void);

/*! Query if mesh library is initialized
    \return true if library is initialized, false if not */
MESH_API bool
mesh_module_is_initialized(void);

/*! Query version of mesh library
    \return Library version */
MESH_API version_t
mesh_module_version(void);

/*! Parse config declarations from JSON buffer
\param buffer Data buffer
\param size Size of data buffer
\param tokens JSON tokens
\param tokens_count Number of JSON tokens */
MESH_API void
mesh_module_parse_config(const char* path, size_t path_size, const char* buffer, size_t size,
                         const struct json_token_t* tokens, size_t tokens_count);

MESH_API mesh_t*
mesh_allocate(size_t expected_vertex_count, size_t expected_triangle_count);

MESH_API void
mesh_initialize(mesh_t* mesh, size_t expected_vertex_count, size_t expected_triangle_count);

MESH_API void
mesh_finalize(mesh_t* mesh);

MESH_API void
mesh_deallocate(mesh_t* mesh);

MESH_API void
mesh_clear(mesh_t* mesh);

MESH_API mesh_t*
mesh_clone(mesh_t* mesh);

MESH_API void
mesh_merge_mesh(mesh_t* mesh, mesh_t* additional);

MESH_API void
mesh_merge_transformed_mesh(mesh_t* mesh, mesh_t* additional, matrix_t transform);

MESH_API void
mesh_set_coordinate_count(mesh_t* mesh, size_t count);

MESH_API void
mesh_reserve_coordinate_count(mesh_t* mesh, size_t count);

MESH_API void
mesh_set_normal_count(mesh_t* mesh, size_t count);

MESH_API void
mesh_reserve_normal_count(mesh_t* mesh, size_t count);

MESH_API void
mesh_set_vertex_count(mesh_t* mesh, size_t count);

MESH_API void
mesh_reserve_vertex_count(mesh_t* mesh, size_t count);

MESH_API void
mesh_set_triangle_count(mesh_t* mesh, size_t count);

MESH_API void
mesh_reserve_triangle_count(mesh_t* mesh, size_t count);

//! Merge coordinates that has same value
MESH_API void
mesh_merge_coordinate(mesh_t* mesh, real merge_coordinate_distance);

// Merge vertices with attributes that has same value. For best result
// merge coordinates first using mesh_merge_coordinate to maximize
// vertex merging
MESH_API void
mesh_merge_vertex(mesh_t* mesh);

// Compact attributes by removing unused attributes and
// remapping indices
MESH_API void
mesh_compact(mesh_t* mesh);

MESH_API void
mesh_calculate_normals(mesh_t* mesh);

MESH_API void
mesh_calculate_tangents(mesh_t* mesh);

MESH_API void
mesh_calculate_bounds(mesh_t* mesh);

MESH_API void
mesh_calculate_topology(mesh_t* mesh);

MESH_API int
mesh_valid(mesh_t* mesh);

MESH_API bool
mesh_triangle_is_degenerate(mesh_t* mesh, mesh_triangle_t* triangle);

MESH_API vector_t
mesh_triangle_normal(mesh_t* mesh, mesh_triangle_t* triangle);

MESH_API void
mesh_dump(mesh_t* mesh, stream_t* stream);

MESH_API mesh_t*
mesh_undump(stream_t* stream);
