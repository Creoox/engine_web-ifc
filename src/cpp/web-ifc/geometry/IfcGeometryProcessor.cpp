/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <spdlog/spdlog.h>

#if defined(DEBUG_DUMP_SVG) || defined(DUMP_CSG_MESHES) || defined(_DEBUG)
#include "../../test/io_helpers.h"
#include "../../test/dumpToThree.h"
#endif

#include "IfcGeometryProcessor.h"
#include <unordered_set>
#include <fstream>
#include <glm/gtx/transform.hpp>
#include "representation/geometry.h"
#include "operations/geometryutils.h"
#include "operations/curve-utils.h"
#include "operations/mesh_utils.h"
#include "operations/meshCleanup.h"
#include "operations/boolean-utils/fuzzy-bools.h"

#if __has_include(<manifold/manifold.h>)
#define WEBIFC_HAS_MANIFOLD 1
#include <manifold/manifold.h>
#else
#define WEBIFC_HAS_MANIFOLD 0
#endif

namespace webifc::geometry
{

    double BOOLSTATUS = 0;

    IfcGeometryProcessor::IfcGeometryProcessor(webifc::parsing::IfcLoader &loader, const webifc::schema::IfcSchemaManager &schemaManager, uint16_t circleSegments, bool coordinateToOrigin, double TOLERANCE_PLANE_INTERSECTION, double TOLERANCE_PLANE_DEVIATION, double TOLERANCE_BACK_DEVIATION_DISTANCE, double TOLERANCE_INSIDE_OUTSIDE_PERIMETER, double TOLERANCE_SCALAR_EQUALITY, double PLANE_REFIT_ITERATIONS)
        : _loader(loader), _cache(loader), _schemaManager(schemaManager), _geometryLoader(loader, _cache, circleSegments)
    {
        _settings._coordinateToOrigin = coordinateToOrigin;
        _settings._circleSegments = circleSegments;
        _settings.TOLERANCE_PLANE_INTERSECTION = TOLERANCE_PLANE_INTERSECTION;
        _settings.TOLERANCE_PLANE_DEVIATION = TOLERANCE_PLANE_DEVIATION;
        _settings.TOLERANCE_BACK_DEVIATION_DISTANCE = TOLERANCE_BACK_DEVIATION_DISTANCE;
        _settings.TOLERANCE_INSIDE_OUTSIDE_PERIMETER = TOLERANCE_INSIDE_OUTSIDE_PERIMETER;
        SetEpsilons(TOLERANCE_SCALAR_EQUALITY, PLANE_REFIT_ITERATIONS);
    }

    IfcGeometryLoader& IfcGeometryProcessor::GetLoader()
    {
         return _geometryLoader;
    }

    void IfcGeometryProcessor::SetTransformation(const std::array<double, 16> &val)
    {
        glm::dmat4 transformation;
        glm::dvec4 v1(val[0], val[1], val[2], val[3]);
        glm::dvec4 v2(val[4], val[5], val[6], val[7]);
        glm::dvec4 v3(val[8], val[9], val[10], val[11]);
        glm::dvec4 v4(val[12], val[13], val[14], val[15]);

        transformation[0] = v1;
        transformation[1] = v2;
        transformation[2] = v3;
        transformation[3] = v4;
        _transformation = transformation;
    }

    IfcGeometry &IfcGeometryProcessor::GetGeometry(uint32_t expressID)
    {
        return _expressIDToGeometry[expressID];
    }

    void IfcGeometryProcessor::Clear()
    {
        std::unordered_map<uint32_t, IfcGeometry>().swap(_expressIDToGeometry);
        _cache.Clear();
    }

    std::array<double, 16> IfcGeometryProcessor::GetFlatCoordinationMatrix() const
    {
        std::array<double, 16> flatTransformation;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                flatTransformation[i * 4 + j] = _coordinationMatrix[i][j];
            }
        }
        return flatTransformation;
    }

    glm::dmat4 IfcGeometryProcessor::GetCoordinationMatrix() const
    {
        return _coordinationMatrix;
    }

    std::optional<glm::dvec4> IfcGeometryProcessor::GetStyleItemFromExpressId(uint32_t expressID)
    {
        std::optional<glm::dvec4> styledItemColor;
        auto &styledItems = _cache.GetStyledItems();
        auto &relMaterials = _cache.GetRelMaterials();
        auto &materialDefinitions = _cache.GetMaterialDefinitions();
        auto styledItem = styledItems.find(expressID);
        if (styledItem != styledItems.end())
        {
            auto items = styledItem->second;
            for (auto item : items)
            {
                styledItemColor = _geometryLoader.GetColor(item.second);
                if (styledItemColor)
                    break;
            }
        }

        if (!styledItemColor)
        {
            auto material = relMaterials.find(expressID);
            if (material != relMaterials.end())
            {
                auto &materials = material->second;
                for (auto item : materials)
                {
                    if (materialDefinitions.count(item.second) != 0)
                    {
                        auto &defs = materialDefinitions.at(item.second);
                        for (auto def : defs)
                        {
                            styledItemColor = _geometryLoader.GetColor(def.second);
                            if (styledItemColor)
                                break;
                        }
                    }

                    // if no color found, check material itself
                    if (!styledItemColor)
                    {
                        styledItemColor = _geometryLoader.GetColor(item.second);
                        if (styledItemColor)
                            break;
                    }
                }
            }
        }
        return styledItemColor;
    }

    IfcComposedMesh IfcGeometryProcessor::GetMesh(uint32_t expressID)
    {
        spdlog::debug("[GetMesh({})]", expressID);
        auto lineType = _loader.GetLineType(expressID);
        auto &relVoids = _cache.GetRelVoids();

        IfcComposedMesh mesh;
        mesh.expressID = expressID;
        std::optional<glm::dvec4> generatedColor = GetStyleItemFromExpressId(expressID);
        if (!generatedColor)
        {
            mesh.color = glm::dvec4(1.0);
            mesh.hasColor = false;
        }
        else
        {
            mesh.color = generatedColor.value();
            mesh.hasColor = true;
        }

        mesh.transformation = glm::dmat4(1);

        if (_schemaManager.IsIfcElement(lineType))
        {
            _loader.MoveToArgumentOffset(expressID, 5);
            uint32_t localPlacement = 0;
            if (_loader.GetTokenType() == parsing::IfcTokenType::REF)
            {
                _loader.StepBack();
                localPlacement = _loader.GetRefArgument();
            }
            uint32_t ifcPresentation = 0;
            if (_loader.GetTokenType() == parsing::IfcTokenType::REF)
            {
                _loader.StepBack();
                ifcPresentation = _loader.GetRefArgument();
            }

            if (localPlacement != 0 && _loader.IsValidExpressID(localPlacement))
            {
                mesh.transformation = _geometryLoader.GetLocalPlacement(localPlacement);
            }

            if (ifcPresentation != 0 && _loader.IsValidExpressID(ifcPresentation))
            {
                mesh.children.push_back(GetMesh(ifcPresentation));
            }

            auto relVoidsIt = relVoids.find(expressID);
            if (relVoidsIt != relVoids.end() && !relVoidsIt->second.empty())
            {
                auto origin = GetOrigin(mesh, _expressIDToGeometry);
                auto normalizeMat = glm::translate(-origin);

                std::vector<IfcGeometry> voidGeoms;

                for (auto relVoidExpressID : relVoidsIt->second)
                {
                    IfcComposedMesh voidGeom = GetMesh(relVoidExpressID);

                    auto flatVoidMesh = flatten(voidGeom, _expressIDToGeometry, normalizeMat);
                    voidGeoms.insert(voidGeoms.end(), flatVoidMesh.begin(), flatVoidMesh.end());
                }

                ApplyBooleanToMeshChildren(mesh, voidGeoms, "DIFFERENCE", _settings, normalizeMat);

#ifdef CSG_DEBUG_OUTPUT
                IfcGeometry testGeo;
                std::vector<IfcGeometry> geomResult = flatten(mesh, _expressIDToGeometry, normalizeMat);
                for (auto &geom : geomResult)
                {
                    testGeo.MergeGeometry(geom);
                }
                webifc::io::DumpIfcGeometry(testGeo, "test.obj");
                std::cout << "Dumped test.obj" << std::endl;
#endif
                return mesh;
            }
            else
            {
                return mesh;
            }

        }
        else
        {
            switch (lineType)
            {
            case schema::IFCSECTIONEDSOLIDHORIZONTAL:
            case schema::IFCSECTIONEDSOLID:
            case schema::IFCSECTIONEDSURFACE:
            {
                auto geom = SectionedSurface(_geometryLoader.GetCrossSections3D(expressID),EPS_SMALL);

                mesh.transformation = glm::dmat4(1);
                
                _expressIDToGeometry[expressID] = geom;
                mesh.hasGeometry = true;

                // #ifdef DEBUG_DUMP_SVG
                //     webifc::io::DumpSectionCurves(curves,"sectioned.obj");
                //     webifc::io::DumpIfcGeometry(geom, "geom.obj");
                // #endif

                return mesh;
            }
            case schema::IFCMAPPEDITEM:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t ifcPresentation = _loader.GetRefArgument();
                uint32_t localPlacement = _loader.GetRefArgument();

                mesh.transformation = _geometryLoader.GetLocalPlacement(localPlacement);
                mesh.children.push_back(GetMesh(ifcPresentation));

                return mesh;
            }
            case schema::IFCBOOLEANCLIPPINGRESULT:
            {
                _loader.MoveToArgumentOffset(expressID, 1);
                uint32_t firstOperandID = _loader.GetRefArgument();
                uint32_t secondOperandID = _loader.GetRefArgument();

                auto firstMesh = GetMesh(firstOperandID);
                auto secondMesh = GetMesh(secondOperandID);

                auto origin = GetOrigin(firstMesh, _expressIDToGeometry);
                auto normalizeMat = glm::translate(-origin);

                auto flatFirstMeshes = flatten(firstMesh, _expressIDToGeometry, normalizeMat);
                auto flatSecondMeshes = flatten(secondMesh, _expressIDToGeometry, normalizeMat);

                IfcGeometry resultMesh = BoolProcess(flatFirstMeshes, flatSecondMeshes, "DIFFERENCE", _settings);

                _expressIDToGeometry[expressID] = resultMesh;
                mesh.hasGeometry = true;
                mesh.transformation = glm::translate(origin);

                if (!mesh.hasColor && firstMesh.hasColor)
                {
                    mesh.hasColor = true;
                    mesh.color = firstMesh.color;
                }

                return mesh;
            }
            case schema::IFCBOOLEANRESULT:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                std::string_view op = _loader.GetStringArgument();

                if (op != "DIFFERENCE" && op != "UNION")
                {
                    spdlog::error("[GetMesh()] Unsupported boolean op {} for entity ID {}", std::string(op), expressID);
                    return mesh;
                }

                uint32_t firstOperandID = _loader.GetRefArgument();
                uint32_t secondOperandID = _loader.GetRefArgument();

                IfcComposedMesh firstMesh = GetMesh(firstOperandID);
                IfcComposedMesh secondMesh = GetMesh(secondOperandID);
                firstMesh.expressID = firstOperandID;
                secondMesh.expressID = secondOperandID;

                auto origin = GetOrigin(firstMesh, _expressIDToGeometry);
                auto normalizeMat = glm::translate(-origin);

                auto flatFirstMeshes = flatten(firstMesh, _expressIDToGeometry, normalizeMat);
                auto flatSecondMeshes = flatten(secondMesh, _expressIDToGeometry, normalizeMat);

                if (flatFirstMeshes.size() == 0)
                {
                    // bail out because we will get strange meshes
                    // if this happens, probably there's an issue parsing the first mesh
                    return mesh;
                }

                size_t numFacesGeoms = 0;
                for (auto& geom : flatFirstMeshes) {
                    numFacesGeoms += geom.numFaces;
                    if (geom.entityID == UINT32_MAX) {
                        geom.entityID = firstOperandID;
                    }
                }
                for (auto& geom : flatSecondMeshes) {
                    numFacesGeoms += geom.numFaces;
                    if (geom.entityID == UINT32_MAX) {
                        geom.entityID = secondOperandID;
                    }
                }

                if (numFacesGeoms > _settings._CSG_MAX_NUM_FACES) {
                    spdlog::warn("High number of faces in CSG operand ({}), skipping further CSG operations", numFacesGeoms);

                    // Merge into a single geometry (equivalent to the first operand without boolean)
                    IfcGeometry resultMesh;
                    for (const auto& g : flatFirstMeshes) {
                        resultMesh.MergeGeometry(g);
                    }
                    resultMesh.entityID = expressID;
                    _expressIDToGeometry[expressID] = resultMesh;
                    mesh.hasGeometry = true;
                    mesh.transformation = glm::translate(origin);

                    if (!mesh.hasColor && firstMesh.hasColor) {
                        mesh.hasColor = true;
                        mesh.color = firstMesh.color;
                    }

                    // No need to set mesh = firstMesh; we're building a flat geometry here for consistency
                    return mesh;
                }

                IfcGeometry resultMesh = BoolProcess(flatFirstMeshes, flatSecondMeshes, std::string(op), _settings);

                resultMesh.entityID = expressID;
                _expressIDToGeometry[expressID] = resultMesh;
                mesh.hasGeometry = true;
                mesh.transformation = glm::translate(origin);
                if (!mesh.hasColor && firstMesh.hasColor)
                {
                    mesh.hasColor = true;
                    mesh.color = firstMesh.color;
                }

                return mesh;
            }
            case schema::IFCHALFSPACESOLID:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t surfaceID = _loader.GetRefArgument();
                std::string_view agreement = _loader.GetStringArgument();

                IfcSurface surface = GetSurface(surfaceID);

                glm::dvec3 extrusionNormal = glm::dvec3(0, 0, 1);

                bool flipWinding = false;
                if (agreement == "T")
                {
                    extrusionNormal *= -1;
                    flipWinding = true;
                }

                double d = EXTRUSION_DISTANCE_HALFSPACE_M / _cache.GetLinearScalingFactor();

                IfcProfile profile;
                profile.isConvex = false;
                profile.curve = GetRectangleCurve(d, d, glm::dmat3(1));

                auto geom = Extrude(profile, extrusionNormal, d);
                geom.halfSpace = true;

                // @Refactor: duplicate of extrudedareasolid
                if (flipWinding)
                {
                    for (uint32_t i = 0; i < geom.numFaces; i++)
                    {
                        uint32_t temp = geom.indexData[i * 3 + 0];
                        geom.indexData[i * 3 + 0] = geom.indexData[i * 3 + 1];
                        geom.indexData[i * 3 + 1] = temp;
                    }
                }

                mesh.transformation = surface.transformation;
                // TODO: this is getting problematic.....
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.hasGeometry = true;

                return mesh;
            }
            case schema::IFCPOLYGONALBOUNDEDHALFSPACE:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t surfaceID = _loader.GetRefArgument();
                std::string_view agreement = _loader.GetStringArgument();
                uint32_t positionID = _loader.GetRefArgument();
                uint32_t boundaryID = _loader.GetRefArgument();

                IfcSurface surface = GetSurface(surfaceID);
                glm::dmat4 position = _geometryLoader.GetLocalPlacement(positionID);
                IfcCurve curve = _geometryLoader.GetCurve(boundaryID, 2);

                if (!curve.IsCCW())
                {
                    curve.Invert();
                }

                glm::dvec3 extrusionNormal = glm::dvec3(0, 0, 1);
                glm::dvec3 planeNormal = surface.transformation[2];
                glm::dvec3 planePosition = surface.transformation[3];

                glm::dmat4 invPosition = glm::inverse(position);
                glm::dvec3 localPlaneNormal = invPosition * glm::dvec4(planeNormal, 0);
                auto localPlanePos = invPosition * glm::dvec4(planePosition, 1);

                bool flipWinding = false;
                double extrudeDistance = EXTRUSION_DISTANCE_HALFSPACE_M / _cache.GetLinearScalingFactor();

                bool halfSpaceInPlaneDirection = agreement != "T";
                bool extrudeInPlaneDirection = glm::dot(localPlaneNormal, extrusionNormal) > 0;
                bool ignoreDistanceInExtrude = (!halfSpaceInPlaneDirection && extrudeInPlaneDirection) || (halfSpaceInPlaneDirection && !extrudeInPlaneDirection);
                if (ignoreDistanceInExtrude)
                {
                    // spec says this should be * 0, but that causes issues for degenerate 0 volume pbhs
                    // hopefully we can get away by just inverting it
                    extrudeDistance *= -1;
                    flipWinding = true;
                }

                IfcProfile profile;
                profile.isConvex = false;
                profile.curve = curve;

                auto geom = Extrude(profile, extrusionNormal, extrudeDistance, localPlaneNormal, localPlanePos);

                if (flipWinding)
                {
                    for (uint32_t i = 0; i < geom.numFaces; i++)
                    {
                        uint32_t temp = geom.indexData[i * 3 + 0];
                        geom.indexData[i * 3 + 0] = geom.indexData[i * 3 + 1];
                        geom.indexData[i * 3 + 1] = temp;
                    }
                }

#ifdef DUMP_CSG_MESHES
                io::DumpIfcGeometry(geom, "pbhs.obj");
#endif

                // TODO: this is getting problematic.....
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.hasGeometry = true;
                mesh.transformation = position;

                return mesh;
            }
            case schema::IFCREPRESENTATIONMAP:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t axis2Placement = _loader.GetRefArgument();
                uint32_t ifcPresentation = _loader.GetRefArgument();

                mesh.transformation = _geometryLoader.GetLocalPlacement(axis2Placement);
                mesh.children.push_back(GetMesh(ifcPresentation));

                return mesh;
            }
            case schema::IFCFACEBASEDSURFACEMODEL:
            case schema::IFCSHELLBASEDSURFACEMODEL:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                auto shells = _loader.GetSetArgument();

                for (auto &shell : shells)
                {
                    uint32_t shellRef = _loader.GetRefArgument(shell);
                    IfcComposedMesh temp;
                    _expressIDToGeometry[shellRef] = GetBrep(shellRef);
                    std::optional<glm::dvec4> shellColor = GetStyleItemFromExpressId(shellRef);
                    if (shellColor)
                    {
                        temp.color = shellColor.value();
                        temp.hasColor = true;
                    }
                    temp.expressID = shellRef;
                    temp.hasGeometry = true;
                    temp.transformation = glm::dmat4(1);
                    mesh.children.push_back(temp);
                }

                int unitaryFaces = 0;
                for (auto &child : mesh.children)
                {
                    auto temp = _expressIDToGeometry[child.expressID];
                    if (temp.numFaces < 4)
                    {
                        unitaryFaces++;
                    }
                }

                IfcGeometry newGeometry;
                newGeometry.entityID = expressID;
                if (unitaryFaces > 12)
                {
                    for (auto &child : mesh.children)
                    {
                        auto temp = _expressIDToGeometry[child.expressID];
                        newGeometry.AddGeometry(temp);
                    }
                    IfcComposedMesh newMesh;
                    _expressIDToGeometry[expressID] = newGeometry;
                    std::optional<glm::dvec4> shellColor = GetStyleItemFromExpressId(expressID);
                    if (shellColor)
                    {
                        newMesh.color = shellColor.value();
                        newMesh.hasColor = true;
                    }
                    newMesh.expressID = expressID;
                    newMesh.hasGeometry = true;
                    newMesh.transformation = glm::dmat4(1);
                    return newMesh;
                }

                return mesh;
            }
            case schema::IFCADVANCEDBREP:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t ifcPresentation = _loader.GetRefArgument();

                _expressIDToGeometry[expressID] = GetBrep(ifcPresentation);
                if (!mesh.hasColor)
                    mesh.color = GetStyleItemFromExpressId(ifcPresentation).value_or(glm::dvec4(1.0));
                mesh.hasGeometry = true;

                return mesh;
            }
            case schema::IFCFACETEDBREP:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t ifcPresentation = _loader.GetRefArgument();

                _expressIDToGeometry[expressID] = GetBrep(ifcPresentation);
                if (!mesh.hasColor)
                    mesh.color = GetStyleItemFromExpressId(ifcPresentation).value_or(glm::dvec4(1.0));
                mesh.hasGeometry = true;

                return mesh;
            }
            case schema::IFCPRODUCTREPRESENTATION:
            case schema::IFCPRODUCTDEFINITIONSHAPE:
            {
                _loader.MoveToArgumentOffset(expressID, 2);
                auto representations = _loader.GetSetArgument();

                for (auto &repToken : representations)
                {
                    uint32_t repID = _loader.GetRefArgument(repToken);
                    mesh.children.push_back(GetMesh(repID));
                }

                return mesh;
            }
            case schema::IFCTOPOLOGYREPRESENTATION:
            case schema::IFCSHAPEREPRESENTATION:
            {
                // IFCTOPOLOGYREPRESENTATION and IFCSHAPEREPRESENTATION are identical in attributes layout
                // attributes: 0=ContextOfItems, 1=RepresentationIdentifier, 2=RepresentationType, 3=Items
                // RepresentationIdentifier is OPTIONAL in the schema, so seek explicitly
                // and token-type check before reading.
                std::string identifierStr;
                _loader.MoveToArgumentOffset(expressID, 1);
                if (_loader.GetTokenType() == parsing::IfcTokenType::STRING) {
                    _loader.StepBack();
                    identifierStr = std::string(_loader.GetStringArgument());
                }

                if (!_settings._representationTypesEnabled.empty()) {
                    if (_settings._representationTypesEnabled.find(identifierStr) == _settings._representationTypesEnabled.end()) {
                        return mesh;
                    }
                }

                _loader.MoveToArgumentOffset(expressID, 3);
                auto repItems = _loader.GetSetArgument();

                for (auto& repToken : repItems)
                {
                    uint32_t repID = _loader.GetRefArgument(repToken);
                    uint32_t repType = _loader.GetLineType(repID);
                    std::string repTypeString = _loader.GetSchemaManager().IfcTypeCodeToType(repType);
                    if (repType == schema::IFCBOUNDINGBOX) {
                        continue;
                    }
                    mesh.children.push_back(GetMesh(repID));
                }

                return mesh;
            }
            case schema::IFCPOLYGONALFACESET:
            {
                _loader.MoveToArgumentOffset(expressID, 0);

                auto coordinatesRef = _loader.GetRefArgument();
                auto points = _geometryLoader.ReadIfcCartesianPointList3D(coordinatesRef);

                // second optional argument closed, ignored

                // indices
                _loader.MoveToArgumentOffset(expressID, 2);
                auto faces = _loader.GetSetArgument();

                IfcGeometry geom;
                geom.entityID = expressID;
                std::vector<IfcBound3D> bounds;
                for (auto &face : faces)
                {
                    uint32_t faceID = _loader.GetRefArgument(face);
                    ReadIndexedPolygonalFace(faceID, bounds, points);

                    TriangulateBounds(geom, bounds, expressID);

                    bounds.clear();
                }

                _loader.MoveToArgumentOffset(expressID, 3);
                if (_loader.GetTokenType() == parsing::IfcTokenType::SET_BEGIN)
                {
                    spdlog::error("[GetMesh()] Unsupported IFCPOLYGONALFACESET with PnIndex {}", expressID);
                }
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.expressID = expressID;
                mesh.hasGeometry = true;

                return mesh;
            }
            case schema::IFCFACESURFACE:
            {
                IfcGeometry geometry;
                _loader.MoveToArgumentOffset(expressID, 0);
                auto bounds = _loader.GetSetArgument();

                std::vector<IfcBound3D> bounds3D(bounds.size());

                for (size_t i = 0; i < bounds.size(); i++)
                {
                    uint32_t boundID = _loader.GetRefArgument(bounds[i]);
                    bounds3D[i] = _geometryLoader.GetBound(boundID);
                }

                _loader.MoveToArgumentOffset(expressID, 1);
                auto surfRef = _loader.GetRefArgument();

                auto surface = GetSurface(surfRef);

                if (surface.BSplineSurface.Active)
                {
                    TriangulateBspline(geometry, bounds3D, surface, _cache.GetLinearScalingFactor());
                }
                else if (surface.CylinderSurface.Active)
                {
                    TriangulateCylindricalSurface(geometry, bounds3D, surface, _settings._circleSegments);
                }
                else if (surface.RevolutionSurface.Active)
                {
                    TriangulateRevolution(geometry, bounds3D, surface, _settings._circleSegments);
                }
                else if (surface.ExtrusionSurface.Active)
                {
                    TriangulateExtrusion(geometry, bounds3D, surface);
                }
                else
                {
                    TriangulateBounds(geometry, bounds3D, expressID);
                }
                geometry.entityID = expressID;
                _expressIDToGeometry[expressID] = geometry;
                mesh.expressID = expressID;
                mesh.hasGeometry = true;

                break;
            }
            case schema::IFCTRIANGULATEDIRREGULARNETWORK:
            case schema::IFCTRIANGULATEDFACESET:
            {
                _loader.MoveToArgumentOffset(expressID, 0);

                auto coordinatesRef = _loader.GetRefArgument();
                auto points = _geometryLoader.ReadIfcCartesianPointList3D(coordinatesRef);

                // second argument normals, ignored
                // third argument closed, ignored

                // indices
                _loader.MoveToArgumentOffset(expressID, 3);
                auto indices = Read2DArrayOfThreeIndices();

                IfcGeometry geom;

                _loader.MoveToArgumentOffset(expressID, 4);
                if (_loader.GetTokenType() == parsing::IfcTokenType::SET_BEGIN)
                {
                    _loader.StepBack();
                    auto pnIndex = Read2DArrayOfThreeIndices();

                    // ignore
                    // std::cout << "Unsupported IFCTRIANGULATEDFACESET with PnIndex!" << std::endl;
                }

                for (size_t i = 0; i < indices.size(); i += 3)
                {
                    int i1 = indices[i + 0] - 1;
                    int i2 = indices[i + 1] - 1;
                    int i3 = indices[i + 2] - 1;

                    geom.AddFace(points[i1], points[i2], points[i3]);
                }

                // DumpIfcGeometry(geom, "test.obj");
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.expressID = expressID;
                mesh.hasGeometry = true;

                return mesh;
            }
            case schema::IFCSURFACECURVESWEPTAREASOLID:
            {

                // TODO: closed sweeps not implemented
                // TODO: the plane is not being used now

                _loader.MoveToArgumentOffset(expressID, 0);

                IfcProfile profile;
                glm::dmat4 placement(1);
                IfcCurve directrix;
                IfcSurface surface;

                double startParam = 0;
                double endParam = 1;
                auto profileID = _loader.GetRefArgument();
                auto placementID = _loader.GetRefArgument();
                auto directrixRef = _loader.GetRefArgument();
                bool closed = false;

                if (_loader.GetTokenType() == parsing::IfcTokenType::REAL)
                {
                    _loader.StepBack();
                    startParam = _loader.GetDoubleArgument();
                }

                if (_loader.GetTokenType() == parsing::IfcTokenType::REAL)
                {
                    _loader.StepBack();
                    endParam = _loader.GetDoubleArgument();
                }

                auto surfaceID = _loader.GetRefArgument();

                if (profileID)
                {
                    profile = _geometryLoader.GetProfile(profileID);
                }
                else
                {
                    break;
                }

                if (placementID)
                {
                    placement = _geometryLoader.GetLocalPlacement(placementID);
                }

                if (directrixRef)
                {
                    directrix = _geometryLoader.GetCurve(directrixRef, 3);
                }
                else
                {
                    break;
                }

                double dst = glm::distance(directrix.points[0], directrix.points[directrix.points.size() - 1]);
                if (startParam == 0 && endParam == 1 && dst < 1e-5)
                {
                    closed = true;
                }

                if (surfaceID)
                {
                    surface = GetSurface(surfaceID);
                }
                else
                {
                    break;
                }

                std::reverse(profile.curve.points.begin(), profile.curve.points.end());

                IfcGeometry geom = Sweep(_cache.GetLinearScalingFactor(), closed, profile, directrix, surface.normal(), true);

                mesh.transformation = placement;
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.expressID = expressID;
                mesh.hasGeometry = true;

                return mesh;
            }
            case schema::IFCFIXEDREFERENCESWEPTAREASOLID:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t profileID = _loader.GetRefArgument();
                uint32_t placementID = _loader.GetOptionalRefArgument();
                uint32_t directrixRef = _loader.GetRefArgument();
                uint32_t fixedReferenceID = _loader.GetRefArgument();

                // Retrieve profile, placement, directrix, and fixed reference direction
                IfcProfile profile = _geometryLoader.GetProfile(profileID);
                glm::dmat4 placement = placementID != UINT32_MAX ? _geometryLoader.GetLocalPlacement(placementID) : glm::dmat4(1.0);
                IfcCurve directrix = _geometryLoader.GetCurve(directrixRef, 3);
                glm::dvec3 fixedReference = _geometryLoader.GetCartesianPoint3D(fixedReferenceID);

                // Check for valid profile and directrix
                if (profile.curve.points.empty() || directrix.points.empty()) {
                    spdlog::error("[GetMesh()] Invalid profile or directrix for IFCFIXEDREFERENCESWEPTAREASOLID {}", expressID);
                    return mesh;
                }

                // Determine if the sweep is closed
                bool closed = glm::distance(directrix.points[0], directrix.points[directrix.points.size() - 1]) < EPS_SMALL;

                // Generate geometry by sweeping the profile with fixed orientation
                IfcGeometry geom = SweepFixedReference(
                    _cache.GetLinearScalingFactor(),
                    closed,
                    profile,
                    directrix,
                    fixedReference
                );

                // Store the geometry and update mesh
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.expressID = expressID;
                mesh.hasGeometry = true;
                mesh.transformation = placement;

                return mesh;
            }
            case schema::IFCSWEPTDISKSOLID:
            {
                // TODO: prevent self intersections in Sweep function still not working properly
                bool closed = false;

                _loader.MoveToArgumentOffset(expressID, 0);
                auto directrixRef = _loader.GetRefArgument();

                double radius = _loader.GetDoubleArgument();
                // double innerRadius = 0.0;

                if (_loader.GetTokenType() == parsing::IfcTokenType::REAL)
                {
                    spdlog::error("[GetMesh()] Inner radius of IFCSWEPTDISKSOLID currently not supported {}", expressID);
                    _loader.StepBack();
                    _loader.GetDoubleArgument();
                }

                // double startParam = 0;
                // double endParam = 0;

                if (_loader.GetTokenType() == parsing::IfcTokenType::REAL)
                {
                    _loader.StepBack();
                    _loader.GetDoubleArgument();
                }

                if (_loader.GetTokenType() == parsing::IfcTokenType::REAL)
                {
                    _loader.StepBack();
                    _loader.GetDoubleArgument();
                }

                IfcCurve directrix = _geometryLoader.GetCurve(directrixRef, 3);

                IfcProfile profile;
                profile.curve = GetCircleCurve(radius, _settings._circleSegments);

                IfcGeometry geom = SweepCircular(_cache.GetLinearScalingFactor(), closed, profile, radius, directrix);

                geom.sweptDiskSolid.axis = std::vector<IfcCurve>{directrix};
                geom.sweptDiskSolid.profiles = std::vector<IfcProfile>{profile};
                geom.sweptDiskSolid.profileRadius = radius;
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.expressID = expressID;
                mesh.hasGeometry = true;

                return mesh;
            }
            case schema::IFCREVOLVEDAREASOLID:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t profileID = _loader.GetRefArgument();
                uint32_t placementID = _loader.GetOptionalRefArgument();
                uint32_t axis1PlacementID = _loader.GetRefArgument();
                double angle = angleConversion(_loader.GetDoubleArgument(), _cache.GetAngleUnits());

                IfcProfile profile = _geometryLoader.GetProfile(profileID);
                glm::dvec3 axisDir = glm::normalize(_geometryLoader.GetAxis1Placement(axis1PlacementID)[0]);
                glm::dvec3 axisPos = _geometryLoader.GetAxis1Placement(axis1PlacementID)[1];

                uint32_t numRots = std::max<uint32_t>(2, _settings._circleSegments);

                // Revolve a closed 2D profile around (axisPos, axisDir) by `angle` radians. The previous implementation built a tiny-radius
                // directrix via BuildArc(pos, ...) where `pos` lies on the axis, so Sweep ran against a degenerate arc and produced
                // missing/pinched geometry. Direct revolution keeps each profile point at its original 3D location at angle 0 and rotates
                // only the perpendicular-to-axis component around the axis. Side surface only; end caps would need a separate step if
                // the revolution is partial and not bounded by other adjacent solids.
                auto revolve = [&](const std::vector<glm::dvec3>& profilePts) -> IfcGeometry {
                    IfcGeometry out;
                    if (profilePts.size() < 2 || numRots < 2) return out;

                    std::vector<std::vector<glm::dvec3>> rings(numRots);
                    for (auto const& p : profilePts) {
                        glm::dvec3 rel = p - axisPos;
                        double along = glm::dot(rel, axisDir);
                        glm::dvec3 parallel = along * axisDir;
                        glm::dvec3 perp = rel - parallel;
                        double r = glm::length(perp);
                        glm::dvec3 u = (r > 1e-12) ? perp / r : glm::dvec3(0);
                        glm::dvec3 v = glm::cross(axisDir, u);
                        for (uint32_t k = 0; k < numRots; ++k) {
                            double t = static_cast<double>(k) / static_cast<double>(numRots - 1);
                            double theta = t * angle;
                            glm::dvec3 pt = axisPos + parallel + r * (std::cos(theta) * u + std::sin(theta) * v);
                            rings[k].push_back(pt);
                        }
                    }

                    for (uint32_t k = 0; k + 1 < numRots; ++k) {
                        auto const& r0 = rings[k];
                        auto const& r1 = rings[k + 1];
                        size_t const n = std::min(r0.size(), r1.size());
                        for (size_t i = 0; i + 1 < n; ++i) {
                            out.AddFace(r0[i], r0[i + 1], r1[i]);
                            out.AddFace(r1[i], r0[i + 1], r1[i + 1]);
                        }
                    }
                    return out;
                };

                IfcGeometry geom;
                if (!profile.isComposite)
                {
                    geom = revolve(profile.curve.points);
                }
                else
                {
                    for (auto const& subProfile : profile.profiles)
                    {
                        IfcGeometry geom_t = revolve(subProfile.curve.points);
                        geom.AddPart(geom_t);
                        geom.AddGeometry(geom_t);
                    }
                }

                if (placementID != UINT32_MAX)
                {
                    mesh.transformation = _geometryLoader.GetLocalPlacement(placementID);
                }
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.expressID = expressID;
                mesh.hasGeometry = true;
                if (!mesh.hasColor)
                    mesh.color = glm::dvec4(1.0);
                else
                    mesh.color = generatedColor.value();
                return mesh;
            }
            case schema::IFCEXTRUDEDAREASOLID:
            case schema::IFCEXTRUDEDAREASOLIDTAPERED:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t profileID = _loader.GetRefArgument();
                uint32_t placementID = _loader.GetOptionalRefArgument();
                uint32_t directionID = _loader.GetRefArgument();
                double depth = _loader.GetDoubleArgument();

                auto lineProfileType = _loader.GetLineType(profileID);
                IfcProfile profile = _geometryLoader.GetProfile(profileID);
                if (!profile.isComposite)
                {
                    if (profile.curve.points.empty())
                    {
                        return mesh;
                    }
                }
                else
                {
                    for (uint32_t i = 0; i < profile.profiles.size(); i++)
                    {
                        if (profile.profiles[i].curve.points.empty())
                        {
                            return mesh;
                        }
                    }
                }

                if (placementID != UINT32_MAX)
                {
                    mesh.transformation = _geometryLoader.GetLocalPlacement(placementID);
                }

                glm::dvec3 dir = _geometryLoader.GetCartesianPoint3D(directionID);

                double dirDot = glm::dot(dir, glm::dvec3(0, 0, 1));
                bool flipWinding = dirDot < 0; // can't be perp according to spec

// TODO: correct dump in case of compositeProfile
#ifdef CSG_DEBUG_OUTPUT
//    io::DumpSVGCurve(profile.curve.points, "IFCEXTRUDEDAREASOLID_curve.html");
#endif

                IfcGeometry geom;

                if (!profile.isComposite)
                {
                    geom = Extrude(profile, dir, depth);
                    if (flipWinding)
                    {
                        for (uint32_t i = 0; i < geom.numFaces; i++)
                        {
                            uint32_t temp = geom.indexData[i * 3 + 0];
                            geom.indexData[i * 3 + 0] = geom.indexData[i * 3 + 1];
                            geom.indexData[i * 3 + 1] = temp;
                        }
                    }
                }
                else
                {
                    for (uint32_t i = 0; i < profile.profiles.size(); i++)
                    {
                        IfcGeometry geom_t = Extrude(profile.profiles[i], dir, depth);
                        if (flipWinding)
                        {
                            for (uint32_t k = 0; k < geom_t.numFaces; k++)
                            {
                                uint32_t temp = geom_t.indexData[k * 3 + 0];
                                geom_t.indexData[k * 3 + 0] = geom_t.indexData[k * 3 + 1];
                                geom_t.indexData[k * 3 + 1] = temp;
                            }
                        }
                        geom.AddPart(geom_t);
                        geom.AddGeometry(geom_t);
                    }
                }

// TODO: correct dump in case of compositeProfile
#ifdef CSG_DEBUG_OUTPUT
//    io::DumpIfcGeometry(geom, "IFCEXTRUDEDAREASOLID_geom.obj");
#endif
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.expressID = expressID;
                mesh.hasGeometry = true;

                return mesh;
            }
            case schema::IFCRIGHTCIRCULARCYLINDER:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t placementID = _loader.GetRefArgument();
                double height = _loader.GetDoubleArgument();
                double radius = _loader.GetDoubleArgument();

                // Create a circular profile
                IfcProfile profile;
                profile.isConvex = true;
                profile.curve = GetCircleCurve(radius, _settings._circleSegments);

                // Extrude along Z-axis
                glm::dvec3 extrusionDir = glm::dvec3(0, 0, 1);
                IfcGeometry geom = Extrude(profile, extrusionDir, height);

                // Set transformation
                if (placementID)
                {
                    mesh.transformation = _geometryLoader.GetLocalPlacement(placementID);
                }

#ifdef CSG_DEBUG_OUTPUT
                io::DumpIfcGeometry(geom, "IFCRIGHTCIRCULARCYLINDER_geom.obj");
#endif
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.expressID = expressID;
                mesh.hasGeometry = true;
                return mesh;
            }
            case schema::IFCGEOMETRICSET:
            case schema::IFCGEOMETRICCURVESET:
            {
                _loader.MoveToArgumentOffset(expressID, 0);
                auto items = _loader.GetSetArgument();

                for (auto &item : items)
                {
                    uint32_t itemID = _loader.GetRefArgument(item);
                    mesh.children.push_back(GetMesh(itemID));
                }

                return mesh;
            }
            case schema::IFCBOUNDINGBOX:
                // ignore bounding box
                return mesh;

            case schema::IFCCARTESIANPOINT:
            {
                // IfcCartesianPoint is derived from IfcRepresentationItem and can be used as representation item directly
                IfcGeometry geom;
                auto point = _geometryLoader.GetCartesianPoint3D(expressID);
                geom.vertexData.push_back(point.x);
                geom.vertexData.push_back(point.y);
                geom.vertexData.push_back(point.z);
                geom.vertexData.push_back(0); // needs to be 6 values per vertex
                geom.vertexData.push_back(0);
                geom.vertexData.push_back(1);
                geom.indexData.push_back(0);

                geom.numPoints = 1;
                geom.primitiveType = fuzzybools::PrimitiveType::POLYLINE;
                mesh.hasGeometry = true;
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;

                return mesh;
            }
            case schema::IFCEDGE:
            case schema::IFCEDGECURVE:
            {
                // IfcEdge is derived from IfcRepresentationItem and can be used as representation item directly
                IfcCurve edge = _geometryLoader.GetEdge(expressID);
                IfcGeometry geom;

                for (uint32_t i = 0; i < edge.points.size(); i++)
                {
                    auto vert = edge.points[i];
                    geom.vertexData.push_back(vert.x);
                    geom.vertexData.push_back(vert.y);
                    geom.vertexData.push_back(vert.z);
                    geom.vertexData.push_back(0); // needs to be 6 values per vertex
                    geom.vertexData.push_back(0);
                    geom.vertexData.push_back(1);
                    geom.indexData.push_back(i);
                }
                geom.numPoints = edge.points.size();
                geom.primitiveType = fuzzybools::PrimitiveType::POLYLINE;
                mesh.hasGeometry = true;
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;

                return mesh;
            }
            case schema::IFCSPHERE:
            {
                // IfcSphere is a CSG solid primitive: center + radius
                // Arguments:
                // 0: Position (IfcAxis2Placement3D) - defines center and local orientation
                // 1: Radius

                _loader.MoveToArgumentOffset(expressID, 0);
                uint32_t placementID = _loader.GetRefArgument();
                double radius = _loader.GetDoubleArgument();

                // Get the placement matrix (center at [3], orientation in columns 0-2)
                glm::dmat4 placement = _geometryLoader.GetLocalPlacement(placementID);

                glm::dvec3 center = glm::dvec3(placement[3]);
                glm::dvec3 xAxis = glm::normalize(glm::dvec3(placement[0]));
                glm::dvec3 yAxis = glm::normalize(glm::dvec3(placement[1]));
                glm::dvec3 zAxis = glm::normalize(glm::dvec3(placement[2]));

                // Tessellate sphere using latitude/longitude grid
                // Use circleSegments for azimuthal (longitude) resolution
                // Use half that for latitudinal resolution (reasonable quality)
                uint32_t azimuthSegments = _settings._circleSegments;
                uint32_t altitudeSegments = _settings._circleSegments / 2;
                if (altitudeSegments < 4) altitudeSegments = 4; // minimum for decent sphere

                IfcGeometry geom;

                // Generate vertices (position + dummy normal for consistency with other cases)
                std::vector<glm::dvec3> vertices;
                vertices.reserve((altitudeSegments + 1) * (azimuthSegments + 1));

                for (uint32_t lat = 0; lat <= altitudeSegments; ++lat)
                {
                    double theta = glm::pi<double>() * lat / altitudeSegments; // 0 (north pole) to pi (south pole)
                    double sinTheta = std::sin(theta);
                    double cosTheta = std::cos(theta);

                    for (uint32_t lon = 0; lon <= azimuthSegments; ++lon)
                    {
                        double phi = 2.0 * glm::pi<double>() * lon / azimuthSegments;
                        double sinPhi = std::sin(phi);
                        double cosPhi = std::cos(phi);

                        // Spherical coordinates to Cartesian (in local space)
                        glm::dvec3 localPos(
                            radius * sinTheta * cosPhi,
                            radius * sinTheta * sinPhi,
                            radius * cosTheta
                        );

                        // Transform to world space using placement axes
                        glm::dvec3 pos = center + localPos.x * xAxis + localPos.y * yAxis + localPos.z * zAxis;

                        vertices.push_back(pos);
                    }
                }

                // Build vertex buffer (6 floats per vertex: pos + dummy normal)
                for (const auto& v : vertices)
                {
                    geom.vertexData.push_back(v.x);
                    geom.vertexData.push_back(v.y);
                    geom.vertexData.push_back(v.z);
                    geom.vertexData.push_back(0.0); // dummy normal X
                    geom.vertexData.push_back(0.0); // dummy normal Y
                    geom.vertexData.push_back(1.0); // dummy normal Z
                }

                // Generate triangle indices (quads -> two triangles, poles handled correctly)
                uint32_t vertsPerRing = azimuthSegments + 1;
                for (uint32_t lat = 0; lat < altitudeSegments; ++lat)
                {
                    uint32_t bottom = lat * vertsPerRing;
                    uint32_t top = bottom + vertsPerRing;

                    for (uint32_t lon = 0; lon < azimuthSegments; ++lon)
                    {
                        uint32_t bl = bottom + lon;     // bottom-left
                        uint32_t br = bottom + lon + 1; // bottom-right
                        uint32_t tl = top + lon;        // top-left
                        uint32_t tr = top + lon + 1;    // top-right

                        // First triangle
                        geom.indexData.push_back(bl);
                        geom.indexData.push_back(br);
                        geom.indexData.push_back(tl);

                        // Second triangle
                        geom.indexData.push_back(br);
                        geom.indexData.push_back(tr);
                        geom.indexData.push_back(tl);
                    }
                }

                geom.numFaces = geom.indexData.size() / 3;
                geom.numPoints = static_cast<uint32_t>(vertices.size());
                geom.primitiveType = fuzzybools::PrimitiveType::TRIANGLES;
                geom.buildPlanes();
                geom.entityID = expressID;
                _expressIDToGeometry[expressID] = geom;
                mesh.hasGeometry = true;
                mesh.expressID = expressID;
                mesh.transformation = placement; // apply the sphere's placement
                
                return mesh;
            }
            case schema::IFCCIRCLE:
            case schema::IFCCOMPOSITECURVE:
            case schema::IFCPOLYLINE:
            case schema::IFCINDEXEDPOLYCURVE:
            case schema::IFCTRIMMEDCURVE:
            case schema::IFCGRADIENTCURVE:
            case schema::IFCCURVESEGMENT:
            {
                auto lineProfileType = _loader.GetLineType(expressID);
                IfcCurve curve = _geometryLoader.GetCurve(expressID, 3, false);

                if (curve.points.size() > 0)
                {
                    IfcGeometry geom;

                    for (uint32_t i = 0; i < curve.points.size(); i++)
                    {
                        auto vert = curve.points[i];
                        geom.vertexData.push_back(vert.x);
                        geom.vertexData.push_back(vert.y);
                        geom.vertexData.push_back(vert.z);
                        geom.vertexData.push_back(0); // needs to be 6 values per vertex
                        geom.vertexData.push_back(0);
                        geom.vertexData.push_back(1);
                        geom.indexData.push_back(i);
                    }
                    geom.numPoints = curve.points.size();
                    geom.primitiveType = fuzzybools::PrimitiveType::POLYLINE;
                    mesh.hasGeometry = true;
                    geom.entityID = expressID;
                    _expressIDToGeometry[expressID] = geom;
                }

                return mesh;
            }
            case schema::IFCTEXTLITERAL:
            case schema::IFCTEXTLITERALWITHEXTENT:
                // TODO: save string of the text literal in IfcComposedMesh
                return mesh;
            default:
                if (lineType != 0) {
                    std::string lineTypeString = _schemaManager.IfcTypeCodeToType(lineType);
                    spdlog::error("[GetMesh()] unexpected mesh type {} for entity ID {}", lineTypeString, expressID);
                }
                break;
            }
        }

        return IfcComposedMesh();
    }

    // Synthesize knots/multiplicities for a plain IFCBSPLINESURFACE (no-knots variant).
    // IFC provides only degree, control points, and a CurveType tag; we derive the
    // knot vector implied by that tag. Without this, Nurbs::init() rejects the surface
    // because knot validation sees empty vectors.
    static void synthesizeBSplineKnots(std::string_view curveType, int degree, size_t numCtrlPts,
                                       std::vector<glm::f64>& knots, std::vector<uint32_t>& mults)
    {
        if (numCtrlPts < static_cast<size_t>(degree + 1)) return;
        size_t const n = numCtrlPts - 1;
        size_t const p = static_cast<size_t>(degree);

        if (curveType == "PIECEWISE_BEZIER_KNOTS" && p > 0 && (n - p) % p == 0)
        {
            size_t const segments = (n - p) / p;
            knots.reserve(segments + 2);
            mults.reserve(segments + 2);
            for (size_t i = 0; i <= segments + 1; ++i) {
                knots.push_back(static_cast<glm::f64>(i));
                mults.push_back((i == 0 || i == segments + 1) ? static_cast<uint32_t>(p + 1)
                                                              : static_cast<uint32_t>(p));
            }
        }
        else if (curveType == "UNIFORM_KNOTS")
        {
            size_t const total = n + p + 2;
            knots.reserve(total);
            mults.reserve(total);
            for (size_t i = 0; i < total; ++i) {
                knots.push_back(static_cast<glm::f64>(i));
                mults.push_back(1);
            }
        }
        else
        {
            // QUASI_UNIFORM_KNOTS, UNSPECIFIED, or unrecognized: clamped uniform.
            size_t const unique = n - p + 2;
            knots.reserve(unique);
            mults.reserve(unique);
            for (size_t i = 0; i < unique; ++i) {
                knots.push_back(static_cast<glm::f64>(i));
                mults.push_back((i == 0 || i == unique - 1) ? static_cast<uint32_t>(p + 1) : 1);
            }
        }
    }

    IfcSurface IfcGeometryProcessor::GetSurface(uint32_t expressID)
    {
        spdlog::debug("[GetSurface({})]", expressID);
        auto lineType = _loader.GetLineType(expressID);

        // TODO: IfcSweptSurface and IfcBSplineSurface still missing
        switch (lineType)
        {
        case schema::IFCPLANE:
        {
            IfcSurface surface;

            _loader.MoveToArgumentOffset(expressID, 0);
            uint32_t locationID = _loader.GetRefArgument();
            surface.transformation = _geometryLoader.GetLocalPlacement(locationID);

            return surface;
        }
        case schema::IFCBSPLINESURFACE:
        {
            IfcSurface surface;

            std::vector<std::vector<glm::vec<3, glm::f64>>> ctrolPts;

            _loader.MoveToArgumentOffset(expressID, 0);
            int Udegree = _loader.GetIntArgument();

            _loader.MoveToArgumentOffset(expressID, 1);
            int Vdegree = _loader.GetIntArgument();

            _loader.MoveToArgumentOffset(expressID, 2);
            auto ctrlPointGroups = _loader.GetSetListArgument();
            for (auto &set : ctrlPointGroups)
            {
                std::vector<glm::vec<3, glm::f64>> list;
                for (auto &token : set)
                {
                    uint32_t pointId = _loader.GetRefArgument(token);
                    list.push_back(_geometryLoader.GetCartesianPoint3D(pointId));
                }
                ctrolPts.push_back(list);
            }

            _loader.MoveToArgumentOffset(expressID, 3);
            auto curveType = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 4);
            auto closedU = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 5);
            auto closedV = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 6);
            auto selfIntersect = _loader.GetStringArgument();

            std::vector<uint32_t> UMultiplicity;
            std::vector<uint32_t> VMultiplicity;
            std::vector<glm::f64> UKnots;
            std::vector<glm::f64> VKnots;
            size_t const numU = ctrolPts.size();
            size_t const numV = ctrolPts.empty() ? 0 : ctrolPts[0].size();
            synthesizeBSplineKnots(curveType, Udegree, numU, UKnots, UMultiplicity);
            synthesizeBSplineKnots(curveType, Vdegree, numV, VKnots, VMultiplicity);

            surface.BSplineSurface.Active = true;
            surface.BSplineSurface.UDegree = Udegree;
            surface.BSplineSurface.VDegree = Vdegree;
            surface.BSplineSurface.ControlPoints = ctrolPts;
            surface.BSplineSurface.ClosedU = closedU;
            surface.BSplineSurface.ClosedV = closedV;
            surface.BSplineSurface.CurveType = curveType;
            surface.BSplineSurface.UMultiplicity = UMultiplicity;
            surface.BSplineSurface.VMultiplicity = VMultiplicity;
            surface.BSplineSurface.UKnots = UKnots;
            surface.BSplineSurface.VKnots = VKnots;

            return surface;
        }
        case schema::IFCBSPLINESURFACEWITHKNOTS:
        {
            IfcSurface surface;

            std::vector<std::vector<glm::vec<3, glm::f64>>> ctrolPts;
            std::vector<uint32_t> UMultiplicity;
            std::vector<uint32_t> VMultiplicity;
            std::vector<glm::f64> UKnots;
            std::vector<glm::f64> VKnots;

            _loader.MoveToArgumentOffset(expressID, 0);
            int Udegree = _loader.GetIntArgument();

            _loader.MoveToArgumentOffset(expressID, 1);
            int Vdegree = _loader.GetIntArgument();

            _loader.MoveToArgumentOffset(expressID, 2);
            auto ctrlPointGroups = _loader.GetSetListArgument();
            for (auto &set : ctrlPointGroups)
            {
                std::vector<glm::vec<3, glm::f64>> list;
                for (auto &token : set)
                {
                    uint32_t pointId = _loader.GetRefArgument(token);
                    list.push_back(_geometryLoader.GetCartesianPoint3D(pointId));
                }
                ctrolPts.push_back(list);
            }

            _loader.MoveToArgumentOffset(expressID, 3);
            auto curveType = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 4);
            auto closedU = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 5);
            auto closedV = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 6);
            auto selfIntersect = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 7);
            auto knotSetU = _loader.GetSetArgument();

            _loader.MoveToArgumentOffset(expressID, 8);
            auto knotSetV = _loader.GetSetArgument();

            _loader.MoveToArgumentOffset(expressID, 9);
            auto indexesSetU = _loader.GetSetArgument();

            _loader.MoveToArgumentOffset(expressID, 10);
            auto indexesSetV = _loader.GetSetArgument();

            for (auto &token : knotSetU)
            {
                UMultiplicity.push_back(_loader.GetIntArgument(token));
            }

            for (auto &token : knotSetV)
            {
                VMultiplicity.push_back(_loader.GetIntArgument(token));
            }

            for (auto &token : indexesSetU)
            {
                UKnots.push_back(_loader.GetDoubleArgument(token));
            }

            for (auto &token : indexesSetV)
            {
                VKnots.push_back(_loader.GetDoubleArgument(token));
            }

            surface.BSplineSurface.Active = true;
            surface.BSplineSurface.UDegree = Udegree;
            surface.BSplineSurface.VDegree = Vdegree;
            surface.BSplineSurface.ControlPoints = ctrolPts;
            surface.BSplineSurface.UMultiplicity = UMultiplicity;
            surface.BSplineSurface.VMultiplicity = VMultiplicity;
            surface.BSplineSurface.UKnots = UKnots;
            surface.BSplineSurface.VKnots = VKnots;

            return surface;

            break;
        }
        case schema::IFCRATIONALBSPLINESURFACEWITHKNOTS:
        {
            IfcSurface surface;

            std::vector<std::vector<glm::vec<3, glm::f64>>> ctrolPts;
            std::vector<std::vector<glm::f64>> weightPts;
            std::vector<uint32_t> UMultiplicity;
            std::vector<uint32_t> VMultiplicity;
            std::vector<glm::f64> UKnots;
            std::vector<glm::f64> VKnots;

            _loader.MoveToArgumentOffset(expressID, 0);
            int Udegree = _loader.GetIntArgument();

            _loader.MoveToArgumentOffset(expressID, 1);
            int Vdegree = _loader.GetIntArgument();

            _loader.MoveToArgumentOffset(expressID, 2);
            auto ctrlPointGroups = _loader.GetSetListArgument();
            for (auto &set : ctrlPointGroups)
            {
                std::vector<glm::vec<3, glm::f64>> list;
                for (auto &token : set)
                {
                    uint32_t pointId = _loader.GetRefArgument(token);
                    list.push_back(_geometryLoader.GetCartesianPoint3D(pointId));
                }
                ctrolPts.push_back(list);
            }

            _loader.MoveToArgumentOffset(expressID, 3);
            auto curveType = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 4);
            auto closedU = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 5);
            auto closedV = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 6);
            auto selfIntersect = _loader.GetStringArgument();

            _loader.MoveToArgumentOffset(expressID, 7);
            auto knotSetU = _loader.GetSetArgument();

            _loader.MoveToArgumentOffset(expressID, 8);
            auto knotSetV = _loader.GetSetArgument();

            _loader.MoveToArgumentOffset(expressID, 9);
            auto indexesSetU = _loader.GetSetArgument();

            _loader.MoveToArgumentOffset(expressID, 10);
            auto indexesSetV = _loader.GetSetArgument();

            _loader.MoveToArgumentOffset(expressID, 12);
            auto weightPointGroups = _loader.GetSetListArgument();
            for (auto &set : weightPointGroups)
            {
                std::vector<glm::f64> list;
                for (auto &token : set)
                {
                    list.push_back(_loader.GetDoubleArgument(token));
                }
                weightPts.push_back(list);
            }

            for (auto &token : knotSetU)
            {
                UMultiplicity.push_back(_loader.GetIntArgument(token));
            }

            for (auto &token : knotSetV)
            {
                VMultiplicity.push_back(_loader.GetIntArgument(token));
            }

            for (auto &token : indexesSetU)
            {
                UKnots.push_back(_loader.GetDoubleArgument(token));
            }

            for (auto &token : indexesSetV)
            {
                VKnots.push_back(_loader.GetDoubleArgument(token));
            }

            // if (UKnots[UKnots.size() - 1] != (int)UKnots[UKnots.size() - 1])
            // {
            //     for (uint32_t i = 0; i < UKnots.size(); i++)
            //     {
            //         UKnots[i] = UKnots[i] * (UKnots.size() - 1) / UKnots[UKnots.size() - 1];
            //     }
            // }

            // if (VKnots[VKnots.size() - 1] != (int)VKnots[VKnots.size() - 1])
            // {
            //     for (uint32_t i = 0; i < VKnots.size(); i++)
            //     {
            //         VKnots[i] = VKnots[i] * (VKnots.size() - 1) / VKnots[VKnots.size() - 1];
            //     }
            // }

            // if (closedU == "T")
            // {
            //  std::vector<std::vector<glm::vec<3, glm::f64>>> newCtrolPts;
            //  for (uint32_t i = 0; i < Udegree; i++)
            //  {
            //      newCtrolPts.push_back(ctrolPts[ctrolPts.size() - 1 + (i - Udegree)]);
            //  }
            //  for (uint32_t s = 0; s < ctrolPts.size(); s++)
            //  {
            //      newCtrolPts.push_back(ctrolPts[s]);
            //  }
            //  ctrolPts = newCtrolPts;
            //  UMultiplicity[0] += Udegree;
            // }

            // if (closedV == "T")
            // {
            //  std::vector<std::vector<glm::vec<3, glm::f64>>> newCtrolPts;
            //  for (uint32_t r = 0; r < ctrolPts.size(); r++)
            //  {
            //      std::vector<glm::vec<3, glm::f64>> newSubList;
            //      for (uint32_t i = 0; i < Vdegree; i++)
            //      {
            //          newSubList.push_back(ctrolPts[r][ctrolPts[r].size() - 1 + (i - Vdegree)]);
            //      }
            //      for (uint32_t s = 0; s < ctrolPts[r].size(); s++)
            //      {
            //          newSubList.push_back(ctrolPts[r][s]);
            //      }
            //      newCtrolPts.push_back(newSubList);
            //  }
            //  ctrolPts = newCtrolPts;
            //  VMultiplicity[0] += Vdegree;
            // }

            surface.BSplineSurface.Active = true;
            surface.BSplineSurface.UDegree = Udegree;
            surface.BSplineSurface.VDegree = Vdegree;
            surface.BSplineSurface.ControlPoints = ctrolPts;
            surface.BSplineSurface.UMultiplicity = UMultiplicity;
            surface.BSplineSurface.VMultiplicity = VMultiplicity;
            surface.BSplineSurface.UKnots = UKnots;
            surface.BSplineSurface.VKnots = VKnots;
            surface.BSplineSurface.Weights = weightPts;

            return surface;

            break;
        }
        case schema::IFCCYLINDRICALSURFACE:
        {
            IfcSurface surface;

            _loader.MoveToArgumentOffset(expressID, 0);
            uint32_t locationID = _loader.GetRefArgument();
            surface.transformation = _geometryLoader.GetLocalPlacement(locationID);

            _loader.MoveToArgumentOffset(expressID, 1);
            double radius = _loader.GetDoubleArgument();

            surface.CylinderSurface.Active = true;
            surface.CylinderSurface.Radius = radius;

            return surface;

            break;
        }
        case schema::IFCSURFACEOFREVOLUTION:
        {
            IfcSurface surface;

            _loader.MoveToArgumentOffset(expressID, 0);
            uint32_t profileID = _loader.GetRefArgument();
            IfcProfile profile = _geometryLoader.GetProfile3D(profileID);

            _loader.MoveToArgumentOffset(expressID, 1);
            if (_loader.GetTokenType() == parsing::IfcTokenType::REF)
            {
                _loader.StepBack();
                uint32_t placementID = _loader.GetRefArgument();
                surface.transformation = _geometryLoader.GetLocalPlacement(placementID);
            }

            _loader.MoveToArgumentOffset(expressID, 2);
            uint32_t locationID = _loader.GetRefArgument();

            surface.RevolutionSurface.Active = true;
            surface.RevolutionSurface.Direction = _geometryLoader.GetLocalPlacement(locationID);
            surface.RevolutionSurface.Profile = profile;

            return surface;

            break;
        }
        case schema::IFCSURFACEOFLINEAREXTRUSION:
        {
            IfcSurface surface;

            _loader.MoveToArgumentOffset(expressID, 0);
            uint32_t profileID = _loader.GetRefArgument();
            IfcProfile profile = _geometryLoader.GetProfile(profileID);

            _loader.MoveToArgumentOffset(expressID, 2);
            uint32_t directionID = _loader.GetRefArgument();
            glm::dvec3 direction = _geometryLoader.GetCartesianPoint3D(directionID);

            _loader.MoveToArgumentOffset(expressID, 3);
            double length = 0;
            if (_loader.GetTokenType() == parsing::IfcTokenType::REAL)
            {
                _loader.StepBack();
                length = _loader.GetDoubleArgument();
            }

            surface.ExtrusionSurface.Active = true;
            surface.ExtrusionSurface.Length = length;
            surface.ExtrusionSurface.Profile = profile;
            surface.ExtrusionSurface.Direction = direction;

            _loader.MoveToArgumentOffset(expressID, 1);
            uint32_t locationID = _loader.GetRefArgument();
            surface.transformation = _geometryLoader.GetLocalPlacement(locationID);

            return surface;

            break;
        }
        default:
            std::string lineTypeString = _schemaManager.IfcTypeCodeToType(lineType);
            spdlog::error("[GetSurface()] unexpected surface type {} for entity ID {}", lineTypeString, expressID);
            break;
        }

        return IfcSurface();
    }

    IfcFlatMesh IfcGeometryProcessor::GetFlatMesh(uint32_t expressID, bool applyLinearScalingFactor)
    {
        spdlog::debug("[GetFlatMesh({})]", expressID);
        IfcFlatMesh flatMesh;
        flatMesh.expressID = expressID;

        IfcComposedMesh composedMesh = GetMesh(expressID);
        if (composedMesh.expressID == UINT32_MAX) {
            composedMesh.expressID = expressID;
        }

        glm::dmat4 mat = glm::dmat4(1);
        if (applyLinearScalingFactor)
        {
            mat = glm::scale(glm::dvec3(_cache.GetLinearScalingFactor()));
        }

        glm::dvec4 color = glm::dvec4(1, 1, 1, 1);
        bool hasColor = false;
        AddComposedMeshToFlatMesh(flatMesh, composedMesh, _transformation * NormalizeIFC * mat, color, hasColor);

        return flatMesh;
    }

    /**
     * Extracts all the geometries and their associated colors and transformations from a composed mesh
     * and adds them to the flat mesh geometries field.
     */
    void IfcGeometryProcessor::AddComposedMeshToFlatMesh(IfcFlatMesh &flatMesh, const IfcComposedMesh &composedMesh, const glm::dmat4 &parentMatrix, const glm::dvec4 &color, bool hasColor)
    {

        glm::dvec4 newParentColor = color;
        bool newHasColor = hasColor;
        glm::dmat4 newMatrix = parentMatrix * composedMesh.transformation;

        if (composedMesh.hasColor && !hasColor)
        {
            newHasColor = true;
            newParentColor = composedMesh.color;
        }

        if (composedMesh.hasGeometry)
        {
            IfcPlacedGeometry geometry;

            if (!_isCoordinated && _settings._coordinateToOrigin)
            {
                auto &geom = _expressIDToGeometry[composedMesh.expressID];
                if (geom.numPoints > 0)
                {
                    auto pt = geom.GetPoint(0);
                    auto transformedPt = newMatrix * glm::dvec4(pt, 1);
                    _coordinationMatrix = glm::translate(-glm::dvec3(transformedPt));
                    _isCoordinated = true;
                }
            }

            auto geom = _expressIDToGeometry[composedMesh.expressID];
            if (geom.isPolygon())
            {
                if (!_settings._exportPolylines)
                {
                    return; // only triangles
                }
            }
            if (geometry.testReverse())
                geom.ReverseFaces();

            if (geom.entityID == UINT32_MAX) {
                geom.entityID = composedMesh.expressID;
            }
            auto translation = glm::dmat4(1.0);

            // #1462 Reports having problems with this line, not sure why this is needed
            translation = geom.Normalize();
            geom.entityID = composedMesh.expressID;
            _expressIDToGeometry[composedMesh.expressID] = geom;

            if (!composedMesh.hasColor)
            {
                geometry.color = newParentColor;
            }
            else
            {
                geometry.color = composedMesh.color;
                newParentColor = composedMesh.color;
                newHasColor = composedMesh.hasColor;
            }

            geometry.transformation = _coordinationMatrix * newMatrix * translation;

            geometry.SetFlatTransformation();
            geometry.geometryExpressID = composedMesh.expressID;

            flatMesh.geometries.push_back(geometry);
        }
        else if (composedMesh.hasColor)
        {
            newParentColor = composedMesh.color;
            newHasColor = composedMesh.hasColor;
        }

        for (auto &c : composedMesh.children)
        {
            AddComposedMeshToFlatMesh(flatMesh, c, newMatrix, newParentColor, newHasColor);
        }
    }


    /**
     * This function traverses an IfcComposedMesh, transforming each child
     * IfcGeometry into world coordinates to perform Boolean operations
     * with the geometries in `secondGroups` (e.g., voids).
     *
     * After performing the Boolean operations in world space, the resulting
     * geometry is transformed back into the local coordinate system of the
     * original IfcComposedMesh and updated in `_expressIDToGeometry`.
     *
     * @param composedMesh The composed mesh whose geometry will be processed.
     * @param secondGroups A vector of IfcGeometry elements to use in Boolean operations (e.g., voids).
     * @param op The Boolean operation to perform ("UNION", "SUBTRACT", "INTERSECT", etc.).
     * @param _settings Settings that control precision, tolerance, or other geometry-specific options.
     * @param mat Transformation matrix representing the parent transformation.
     */
    void IfcGeometryProcessor::ApplyBooleanToMeshChildren(IfcComposedMesh &composedMesh, std::vector<IfcGeometry> &secondGroups, std::string op, IfcGeometrySettings _settings, glm::dmat4 mat = glm::dmat4(1))
    {
        const glm::dmat4 newMat = mat * composedMesh.transformation;
        const glm::dmat4 invMat = glm::inverse(newMat);
        bool transformationBreaksWinding = MatrixFlipsTriangles(newMat);
        bool transformationBreaksWindingInverse = MatrixFlipsTriangles(invMat);
        auto geomIt = _expressIDToGeometry.find(composedMesh.expressID);
        if (composedMesh.hasGeometry && geomIt != _expressIDToGeometry.end())
        {
            IfcGeometry meshGeom = geomIt->second;
#ifdef _DEBUG_PRINT
            std::cout << "[ApplyBoolean] expressID=" << composedMesh.expressID << " faces=" << meshGeom.numFaces << " children=" << composedMesh.children.size() << std::endl;
#endif
            if (meshGeom.numFaces <= _settings._CSG_MAX_NUM_FACES)
            {
                std::vector<IfcGeometry> transformedGeoms = transformIfcGeometry(meshGeom, newMat, transformationBreaksWinding);
                IfcGeometry geometryResult = BoolProcess(transformedGeoms, secondGroups, op, _settings);
                
                std::vector<IfcGeometry> localGeom = transformIfcGeometry(geometryResult, invMat, transformationBreaksWindingInverse);
                IfcGeometry localGeomMerged;
                for (const auto& geom : localGeom)
                {
                    localGeomMerged.MergeGeometry(geom);
                }
                _expressIDToGeometry[composedMesh.expressID] = localGeomMerged;
            }
        }
        for (auto &c : composedMesh.children)
        {
            ApplyBooleanToMeshChildren(c, secondGroups, op, _settings, newMat);
        }
    }

    std::vector<uint32_t> IfcGeometryProcessor::Read2DArrayOfThreeIndices()
    {
        std::vector<uint32_t> result;

        _loader.GetTokenType();

        // while we have point set begin
        while (_loader.GetTokenType() == parsing::IfcTokenType::SET_BEGIN)
        {
            result.push_back((uint32_t)_loader.GetIntArgument());
            result.push_back((uint32_t)_loader.GetIntArgument());
            result.push_back((uint32_t)_loader.GetIntArgument());

            // read point set end
            _loader.GetTokenType();
        }

        return result;
    }

    void IfcGeometryProcessor::ReadIndexedPolygonalFace(uint32_t expressID, std::vector<IfcBound3D> &bounds, const std::vector<glm::dvec3> &points)
    {

        spdlog::debug("[ReadIndexedPolygonalFace({})]", expressID);
        auto lineType = _loader.GetLineType(expressID);

        bounds.emplace_back();

        switch (lineType)
        {
        case schema::IFCINDEXEDPOLYGONALFACEWITHVOIDS:
        case schema::IFCINDEXEDPOLYGONALFACE:
        {
            _loader.MoveToArgumentOffset(expressID, 0);
            auto indexIDs = _loader.GetSetArgument();

            IfcGeometry geometry;
            for (auto &indexID : indexIDs)
            {
                uint32_t index = _loader.GetIntArgument(indexID);
                glm::dvec3 point = points[index - 1]; // indices are 1-based

                // I am not proud of this
                bounds.back().curve.points.push_back(point);
            }

            if (lineType == schema::IFCINDEXEDPOLYGONALFACE)
            {
                break;
            }

            // case IFCINDEXEDPOLYGONALFACEWITHVOIDS
            _loader.MoveToArgumentOffset(expressID, 1);

            // guaranteed to be set begin
            _loader.GetTokenType();

            // while we have hole-index set begin
            while (_loader.GetTokenType() == parsing::IfcTokenType::SET_BEGIN)
            {
                bounds.emplace_back();

                while (_loader.GetTokenType() != parsing::IfcTokenType::SET_END)
                {
                    _loader.StepBack();
                    uint32_t index = _loader.GetIntArgument();

                    glm::dvec3 point = points[index - 1]; // indices are still 1-based

                    // I am also not proud of this
                    bounds.back().curve.points.push_back(point);
                }
            }

            break;
        }
        default:
            spdlog::error("[ReadIndexedPolygonalFace()] unexpected indexedface type {}: {}", expressID, lineType);
            break;
        }
    }

    IfcGeometry IfcGeometryProcessor::GetBrep(uint32_t expressID)
    {
        spdlog::debug("[GetBrep({})]", expressID);
        auto lineType = _loader.GetLineType(expressID);
        switch (lineType)
        {
        case schema::IFCCONNECTEDFACESET:
        case schema::IFCCLOSEDSHELL:
        case schema::IFCOPENSHELL:
        {
            _loader.MoveToArgumentOffset(expressID, 0);
            auto faces = _loader.GetSetArgument();

            IfcGeometry geometry;
            for (auto &faceToken : faces)
            {
                uint32_t faceID = _loader.GetRefArgument(faceToken);
                AddFaceToGeometry(faceID, geometry);
            }
            geometry.entityID = expressID;
#ifdef DUMP_CSG_MESHES
            if( lineType == schema::IFCCLOSEDSHELL )
            {
                fuzzybools::Geometry fuzzyGeom = geometry;
                auto meshInfo = meshCleanup::isMeshWatertight(fuzzyGeom);
                if (!meshInfo.watertight) {
                    webifc::io::DumpIfcGeometry(geometry, "IFCCLOSEDSHELL-notWaterTight.obj");
                    spdlog::warn("[GetBrep()] geometry with expressID {} is not watertight after triangulation", expressID);
                }
			}
            
#endif
            return geometry;
        }
        default:
            spdlog::error("[GetBrep()] unexpected shell type {}: {}", expressID, lineType);
            break;
        }

        return IfcGeometry();
    }

    void IfcGeometryProcessor::AddFaceToGeometry(uint32_t expressID, IfcGeometry &geometry)
    {
        spdlog::debug("[AddFaceToGeometry({})]", expressID);
        auto lineType = _loader.GetLineType(expressID);

        switch (lineType)
        {
        case schema::IFCFACE:
        {
            _loader.MoveToArgumentOffset(expressID, 0);
            auto bounds = _loader.GetSetArgument();

            std::vector<IfcBound3D> bounds3D(bounds.size());

            for (size_t i = 0; i < bounds.size(); i++)
            {
                uint32_t boundID = _loader.GetRefArgument(bounds[i]);
                bounds3D[i] = _geometryLoader.GetBound(boundID);
            }

            TriangulateBounds(geometry, bounds3D, expressID);
            break;
        }
        case schema::IFCADVANCEDFACE:
        {
            _loader.MoveToArgumentOffset(expressID, 0);
            auto bounds = _loader.GetSetArgument();

            std::vector<IfcBound3D> bounds3D(bounds.size());

            for (size_t i = 0; i < bounds.size(); i++)
            {
                uint32_t boundID = _loader.GetRefArgument(bounds[i]);
                bounds3D[i] = _geometryLoader.GetBound(boundID);
            }

            _loader.MoveToArgumentOffset(expressID, 1);
            auto surfRef = _loader.GetRefArgument();

            auto surface = GetSurface(surfRef);

            // Read SameSense flag (arg 2): .T. means face normal agrees with surface parametric
            // normal; .F. means the face normal is opposite – winding must be reversed.
            _loader.MoveToArgumentOffset(expressID, 2);
            std::string_view senseStr = _loader.GetStringArgument();
            bool sameSense = (senseStr == "T");

            uint32_t facesBefore = geometry.numFaces;

            // TODO: place the face in the surface and tringulate

            if (surface.BSplineSurface.Active)
            {
                TriangulateBspline(geometry, bounds3D, surface, _cache.GetLinearScalingFactor());
            }
            else if (surface.CylinderSurface.Active)
            {
                TriangulateCylindricalSurface(geometry, bounds3D, surface, _settings._circleSegments);
            }
            else if (surface.RevolutionSurface.Active)
            {
                TriangulateRevolution(geometry, bounds3D, surface, _settings._circleSegments);
            }
            else if (surface.ExtrusionSurface.Active)
            {
                TriangulateExtrusion(geometry, bounds3D, surface);
            }
            else
            {
                TriangulateBounds(geometry, bounds3D, expressID);
            }

            // If SameSense is .F., the face outward normal is opposite to the surface parametric
            // normal, so reverse the winding of all faces added by the triangulation above.
            if (!sameSense)
            {
                for (uint32_t fi = facesBefore; fi < geometry.numFaces; fi++)
                {
                    std::swap(geometry.indexData[fi * 3 + 0], geometry.indexData[fi * 3 + 2]);
                }
            }
            break;
        }
        default:
            spdlog::error("[AddFaceToGeometry()] unexpected face type {}: {}", expressID, lineType);
            break;
        }
    }

    IfcGeometry IfcGeometryProcessor::BoolProcess(const std::vector<IfcGeometry> &firstGeoms, std::vector<IfcGeometry> &secondGeoms, std::string op, IfcGeometrySettings _settings)
    {
        spdlog::debug("[BoolProcess({})]");
        IfcGeometry finalResult;

        for (auto &firstGeom : firstGeoms)
        {
            IfcGeometry firstOperand = firstGeom;
            for (auto &secondGeom : secondGeoms)
            {
                if (secondGeom.numFaces == 0)
                {
                    spdlog::error("[BoolProcess()] bool aborted due to empty source or target");

                    // bail out because we will get strange meshes
                    // if this happens, probably there's an issue parsing the mesh that occurred earlier
                    continue;
                }

                if (firstOperand.numFaces == 0 && op != "UNION")
                {
                    spdlog::error("[BoolProcess()] bool aborted due to empty source or target");

                    // bail out because we will get strange meshes
                    // if this happens, probably there's an issue parsing the mesh that occurred earlier
                    break;
                }

                IfcGeometry secondOperand;

                if (secondGeom.halfSpace)
                {
                    glm::dvec3 origin = secondGeom.halfSpaceOrigin;
                    glm::dvec3 x = secondGeom.halfSpaceX - origin;
                    glm::dvec3 y = secondGeom.halfSpaceY - origin;
                    glm::dvec3 z = secondGeom.halfSpaceZ - origin;
                    glm::dmat4 trans = glm::dmat4(
                        glm::dvec4(x, 0),
                        glm::dvec4(y, 0),
                        glm::dvec4(z, 0),
                        glm::dvec4(0, 0, 0, 1));

                    double scaleX = 1;
                    double scaleY = 1;
                    double scaleZ = 1;

                    for (uint32_t i = 0; i < firstOperand.numPoints; i++)
                    {
                        glm::dvec3 p = firstOperand.GetPoint(i);
                        glm::dvec3 vec = (p - origin);
                        double dx = glm::dot(vec, x);
                        double dy = glm::dot(vec, y);
                        double dz = glm::dot(vec, z);
                        if (glm::abs(dx) > scaleX)
                        {
                            scaleX = glm::abs(dx);
                        }
                        if (glm::abs(dy) > scaleY)
                        {
                            scaleY = glm::abs(dy);
                        }
                        if (glm::abs(dz) > scaleZ)
                        {
                            scaleZ = glm::abs(dz);
                        }
                    }
                    secondOperand.AddGeometry(secondGeom, trans, scaleX * 2, scaleY * 2, scaleZ * 2, secondGeom.halfSpaceOrigin);
                }
                else
                {
                    secondOperand = secondGeom;
                }

                // Pre-boolean operand scale-up: considered and not applied.
                // Empirically the scale-up trick trades one class of artifact for another (window-opening membranes 
                // disappear but edge-bleed membranes appear), so we leave the operand untouched and clean up post hoc.

                firstOperand.buildPlanes();
                secondOperand.buildPlanes();

                fuzzybools::SetEpsilons(_settings.TOLERANCE_PLANE_INTERSECTION, _settings.TOLERANCE_PLANE_DEVIATION, _settings.TOLERANCE_BACK_DEVIATION_DISTANCE, _settings.TOLERANCE_INSIDE_OUTSIDE_PERIMETER, _settings.TOLERANCE_BOUNDING_BOX, BOOLSTATUS);

                if (op == "DIFFERENCE")
                {
                    firstOperand = Subtract(firstOperand, secondOperand);
                }
                else if (op == "UNION")
                {
                    firstOperand = Union(firstOperand, secondOperand);
                }
#ifdef _DEBUG_PRINT
                std::cout << "[BoolProcess] result.faces=" << firstOperand.numFaces << std::endl;
#endif
            }
            finalResult.AddGeometry(firstOperand);
        }

        return finalResult;
    }

    // Safety guard for pathological single-op growth: if a boolean returns a face count that is more 
    // than `kPathologicalFactor` times the sum of its input face counts, something went badly wrong in the kernel.
    // Run a cheap collapse pass (degenerate + exact duplicate triangles) before the result is allowed into the 
    // next op in a chain. This caps compounding damage when the kernel misclassifies in one op.
    static constexpr uint32_t kPathologicalFactor = 10;

    static void CheapCollapseAfterBoolean(fuzzybools::Geometry& g)
    {
        const uint32_t nF = g.numFaces;
        if (nF == 0) return;
        // Exact triangle deduplication using the fvertexData / vertexData positions (no tolerance; only faces that are bit-for-bit 
        // equal after the boolean are dropped).
        struct TriKey { uint64_t a, b, c; };
        struct TriKeyHash {
            size_t operator()(const TriKey& k) const {
                size_t h = std::hash<uint64_t>()(k.a);
                h ^= std::hash<uint64_t>()(k.b) + 0x9e3779b9 + (h << 6) + (h >> 2);
                h ^= std::hash<uint64_t>()(k.c) + 0x9e3779b9 + (h << 6) + (h >> 2);
                return h;
            }
        };
        struct TriKeyEq {
            bool operator()(const TriKey& a, const TriKey& b) const {
                return a.a == b.a && a.b == b.b && a.c == b.c;
            }
        };
        std::unordered_set<TriKey, TriKeyHash, TriKeyEq> seen;
        seen.reserve(nF);
        std::vector<uint32_t> newIdx;
        std::vector<uint32_t> newPlane;
        newIdx.reserve(g.indexData.size());
        newPlane.reserve(g.planeData.size());
        for (uint32_t i = 0; i < nF; i++) {
            uint32_t i0 = g.indexData[i * 3 + 0];
            uint32_t i1 = g.indexData[i * 3 + 1];
            uint32_t i2 = g.indexData[i * 3 + 2];
            if (i0 == i1 || i1 == i2 || i0 == i2) continue; // degenerate
            // Canonicalise the triple so cyclic rotations match.
            uint64_t a = i0, b = i1, c = i2;
            if (b < a && b < c) { uint64_t t = a; a = b; b = c; c = t; }
            else if (c < a && c < b) { uint64_t t = a; a = c; c = b; b = t; }
            if (!seen.insert({ a, b, c }).second) continue; // exact duplicate
            newIdx.push_back(i0);
            newIdx.push_back(i1);
            newIdx.push_back(i2);
            if (i < g.planeData.size()) newPlane.push_back(g.planeData[i]);
        }
        g.indexData = std::move(newIdx);
        g.planeData = std::move(newPlane);
        g.numFaces = (uint32_t)(g.indexData.size() / 3);
        g.hasPlanes = false;
    }

    IfcGeometry IfcGeometryProcessor::Union(IfcGeometry firstOperand, IfcGeometry secondOperand)
    {
        const uint32_t inputFaces = firstOperand.numFaces + secondOperand.numFaces;
        fuzzybools::Geometry result = fuzzybools::Union(firstOperand, secondOperand);
        if (inputFaces > 0 && result.numFaces > kPathologicalFactor * inputFaces)
        {
            CheapCollapseAfterBoolean(result);
        }
        meshCleanup::PostBooleanOperationMeshCleanup(result);
        IfcGeometry out;
        static_cast<fuzzybools::Geometry&>(out) = std::move(result);
        return out;
    }

#if WEBIFC_HAS_MANIFOLD
    // Snap a metre-scale coordinate to a fixed sub-nanometre grid so that ULP-scale numerical noise
    // (~1e-15..1e-16 at metre scale) collapses to exact coincidence. Real CAD features are mm-scale
    // or coarser, so a 1 nm grid is far below anything geometrically meaningful in IFC. Without this,
    // an opening whose extrusion overshoots a wall by 1 ULP produces a 1 ULP-thick residual membrane
    // in manifold's Boolean output -- manifold is exact and treats the gap as real geometry.
    static constexpr double kManifoldSnapGrid = 1.0e-9;
    static inline double SnapForManifold(double v) {
        return std::round(v / kManifoldSnapGrid) * kManifoldSnapGrid;
    }

    static manifold::MeshGL64 ToManifoldMesh(const fuzzybools::Geometry& g)
    {
        manifold::MeshGL64 mesh;
        mesh.numProp = 3;
        mesh.vertProperties.reserve(static_cast<size_t>(g.numPoints) * 3);
        for (uint32_t i = 0; i < g.numPoints; ++i) {
            glm::dvec3 p = g.GetPoint(i);
            mesh.vertProperties.push_back(SnapForManifold(p.x));
            mesh.vertProperties.push_back(SnapForManifold(p.y));
            mesh.vertProperties.push_back(SnapForManifold(p.z));
        }
        mesh.triVerts.reserve(g.indexData.size());
        for (uint32_t idx : g.indexData) mesh.triVerts.push_back(static_cast<uint64_t>(idx));

        // Manifold requires CCW winding viewed from outside (positive signed volume for a closed solid).
        // fuzzybools doesn't enforce this convention, so detect inverted input via the divergence-theorem
        // volume integral and flip per-triangle winding before handoff. Subtract a reference point to keep
        // the sum origin-independent.
        if (g.numFaces > 0) {
            glm::dvec3 ref = g.GetPoint(g.indexData[0]);
            double signedVol = 0.0;
            for (uint32_t i = 0; i < g.numFaces; ++i) {
                glm::dvec3 a = g.GetPoint(g.indexData[3 * i + 0]) - ref;
                glm::dvec3 b = g.GetPoint(g.indexData[3 * i + 1]) - ref;
                glm::dvec3 c = g.GetPoint(g.indexData[3 * i + 2]) - ref;
                signedVol += glm::dot(a, glm::cross(b, c)) / 6.0;
            }
            if (signedVol < 0.0) {
                spdlog::debug("[ToManifoldMesh] numFaces={} signedVol={} (negative => CW input, flipping)", g.numFaces, signedVol);
                for (size_t t = 0; t < mesh.triVerts.size(); t += 3) {
                    std::swap(mesh.triVerts[t + 1], mesh.triVerts[t + 2]);
                }
            }
        }

        mesh.Merge();
        return mesh;
    }

    static fuzzybools::Geometry FromManifold(const manifold::Manifold& m)
    {
        fuzzybools::Geometry out;
        if (m.IsEmpty()) return out;
        manifold::MeshGL64 mesh = m.GetMeshGL64();
        const size_t nTri = mesh.NumTri();
        for (size_t t = 0; t < nTri; ++t) {
            auto v0 = mesh.GetVertPos(mesh.triVerts[3 * t + 0]);
            auto v1 = mesh.GetVertPos(mesh.triVerts[3 * t + 1]);
            auto v2 = mesh.GetVertPos(mesh.triVerts[3 * t + 2]);
            out.AddFace(glm::dvec3(v0.x, v0.y, v0.z),
                        glm::dvec3(v1.x, v1.y, v1.z),
                        glm::dvec3(v2.x, v2.y, v2.z));
        }
        return out;
    }

    // Diagnostic: dump a manifold MeshGL64 as a wavefront OBJ. Positions only, 1-indexed faces, full
    // double precision so coincident planes are byte-comparable in the file. Intended for one-shot
    // CSG investigations; remove or gate behind a flag once diagnosis is complete.
    static void DumpManifoldMeshToObj(const manifold::MeshGL64& mesh, const std::string& path)
    {
        std::ofstream out(path);
        if (!out.is_open()) {
            spdlog::warn("[DumpManifoldMeshToObj] failed to open {}", path);
            return;
        }
        out << std::scientific;
        out.precision(17);
        const size_t nVert = static_cast<size_t>(mesh.NumVert());
        const size_t nTri = static_cast<size_t>(mesh.NumTri());
        for (size_t v = 0; v < nVert; ++v) {
            const size_t off = v * mesh.numProp;
            out << "v " << mesh.vertProperties[off + 0]
                << " "  << mesh.vertProperties[off + 1]
                << " "  << mesh.vertProperties[off + 2] << "\n";
        }
        for (size_t t = 0; t < nTri; ++t) {
            out << "f " << (mesh.triVerts[3 * t + 0] + 1)
                << " "  << (mesh.triVerts[3 * t + 1] + 1)
                << " "  << (mesh.triVerts[3 * t + 2] + 1) << "\n";
        }
    }
#endif

    IfcGeometry IfcGeometryProcessor::Subtract(IfcGeometry firstOperand, IfcGeometry secondOperand)
    {
        const uint32_t inputFaces = firstOperand.numFaces + secondOperand.numFaces;
        fuzzybools::Geometry result;
        bool gotResult = false;

#if WEBIFC_HAS_MANIFOLD
        try {
            manifold::Manifold a(ToManifoldMesh(firstOperand));
            manifold::Manifold b(ToManifoldMesh(secondOperand));
            if (a.Status() == manifold::Manifold::Error::NoError && b.Status() == manifold::Manifold::Error::NoError) {
                manifold::Manifold diff = a.Boolean(b, manifold::OpType::Subtract);
                if (diff.Status() == manifold::Manifold::Error::NoError) {
#ifdef _DEBUG
                    static int kManifoldDumpCounter = 0;
                    const int n = kManifoldDumpCounter++;
                    const std::string base = "manifold_dump_" + std::to_string(n);
                    DumpManifoldMeshToObj(a.GetMeshGL64(),    base + "_a.obj");
                    DumpManifoldMeshToObj(b.GetMeshGL64(),    base + "_b.obj");
                    DumpManifoldMeshToObj(diff.GetMeshGL64(), base + "_diff.obj");
                    spdlog::info("[Subtract] dumped manifold meshes to {}_(a|b|diff).obj (cwd)", base);
#endif
                    result = FromManifold(diff);
                    result.mBoolOpCount = std::max(firstOperand.mBoolOpCount, secondOperand.mBoolOpCount) + 1;
                    gotResult = true;
                }
            }
        } catch (...) {
            gotResult = false;
        }
#endif

        if (!gotResult) {
            result = fuzzybools::Subtract(firstOperand, secondOperand);
        }

        if (inputFaces > 0 && result.numFaces > kPathologicalFactor * inputFaces)
        {
            CheapCollapseAfterBoolean(result);
        }
        meshCleanup::PostBooleanOperationMeshCleanup(result);
        IfcGeometry out;
        static_cast<fuzzybools::Geometry&>(out) = std::move(result);
        return out;
    }

    IfcGeometryProcessor *IfcGeometryProcessor::Clone(const webifc::parsing::IfcLoader &newLoader) const
    {
        IfcGeometryProcessor *newProcessor = new IfcGeometryProcessor(_settings, _expressIDToGeometry, _transformation, newLoader, _schemaManager, _isCoordinated, _expressIdCyl, _expressIdRect, _coordinationMatrix, _predefinedCylinder, _predefinedCube);
        return newProcessor;
    }

    IfcGeometryProcessor::IfcGeometryProcessor(const IfcGeometrySettings &settings, std::unordered_map<uint32_t, IfcGeometry> expressIDToGeometry, glm::dmat4 transformation, const parsing::IfcLoader &loader, const schema::IfcSchemaManager &schemaManager, bool isCoordinated, uint32_t expressIdCyl, uint32_t expressIdRect, glm::dmat4 coordinationMatrix, IfcGeometry predefinedCylinder, IfcGeometry predefinedCube)
        : _settings(settings), _expressIDToGeometry(expressIDToGeometry), _transformation(transformation), _loader(loader), _cache(cache::IfcCache(loader)), _schemaManager(schemaManager), _isCoordinated(isCoordinated), _expressIdCyl(expressIdCyl), _expressIdRect(expressIdRect), _coordinationMatrix(coordinationMatrix), _predefinedCylinder(predefinedCylinder), _predefinedCube(predefinedCube), _geometryLoader(geometry::IfcGeometryLoader(loader,_cache,settings._circleSegments))
    {
    }

}
