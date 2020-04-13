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
