#include <core/core.hpp>
#include <physics/physics_statics.hpp>
#include <rendering/debugrendering.hpp>

namespace legion::physics
{
    bool PhysicsStatics::FindSeperatingAxisByExtremePointProjection(ConvexCollider* convexA
        , ConvexCollider* convexB, const math::mat4& transformA, const math::mat4& transformB,
        PointerEncapsulator<HalfEdgeFace>& refFace,
        float& maximumSeperation, bool shouldDebug )
    {
        float currentMaximumSeperation = std::numeric_limits<float>::lowest();

        for (auto face : convexB->GetHalfEdgeFaces())
        {

            math::vec3 seperatingAxis = math::normalize(transformB * math::vec4((face->normal), 0));

            math::vec3 transformedPositionB = transformB * math::vec4(face->centroid, 1);

            //get extreme point of other face in normal direction
            math::vec3 worldSupportPoint; 
            GetSupportPointNoTransform(transformedPositionB, -seperatingAxis,
                convexA, transformA, worldSupportPoint);

            float seperation = math::dot(worldSupportPoint - transformedPositionB, seperatingAxis);

            if (seperation > currentMaximumSeperation)
            {
                currentMaximumSeperation = seperation;
                refFace.ptr = face;
            }

            if (seperation > 0)
            {
                //we have found a seperating axis, we can exit early
                maximumSeperation = currentMaximumSeperation;
                return true;
            }
        }

        //no seperating axis was found
        maximumSeperation = currentMaximumSeperation;

        return false;
    }

    float PhysicsStatics::GetSupportPoint(const std::vector<math::vec3>& vertices, const math::vec3& direction,math::vec3& outVec)
    {
        float currentMaximumSupportPoint = std::numeric_limits<float>::lowest();

        for (const auto& vert : vertices)
        {
            float dotResult = math::dot(direction, vert);

            if (dotResult > currentMaximumSupportPoint)
            {
                currentMaximumSupportPoint = dotResult;
                outVec = vert;
            }
        }

        return currentMaximumSupportPoint;
    }

    void PhysicsStatics::GetSupportPointNoTransform(math::vec3 planePosition, math::vec3 direction, ConvexCollider* collider, const math::mat4& colliderTransform
        , math::vec3& worldSupportPoint)
    {
        float largestDistanceInDirection = std::numeric_limits<float>::lowest();
        planePosition = math::inverse(colliderTransform) * math::vec4(planePosition, 1);
        direction = math::inverse(colliderTransform) * math::vec4(direction, 0);

        const auto& vertices = collider->GetVertices();

        for (int vertexIndex = 0 ; vertexIndex < vertices.size(); ++vertexIndex)
        {
            math::vec3 transformedVert = math::vec4(vertices.at(vertexIndex), 1);

            float dotResult = math::dot(transformedVert - planePosition, direction);

            if (dotResult > largestDistanceInDirection)
            {
                largestDistanceInDirection = dotResult;
                worldSupportPoint = transformedVert;
            }
        }

        worldSupportPoint = colliderTransform * math::vec4(worldSupportPoint, 1);
    }

    bool PhysicsStatics::FindSeperatingAxisByGaussMapEdgeCheck(ConvexCollider* convexA, ConvexCollider* convexB, const math::mat4& transformA,
        const math::mat4& transformB, PointerEncapsulator<HalfEdgeEdge>& refEdge, PointerEncapsulator<HalfEdgeEdge>& incEdge,
        math::vec3& seperatingAxisFound, float& maximumSeperation, bool shouldDebug)
    {
        float currentMaximumSeperation = std::numeric_limits<float>::lowest();

        math::vec3 centroidDir = transformA * math::vec4(convexA->GetLocalCentroid(), 0);
        math::vec3 positionA = math::vec3(transformA[3]) + centroidDir;

        int facei = 0;
        int facej = 0;

        for (const auto faceA : convexA->GetHalfEdgeFaces())
        {
            //----------------- Get all edges of faceA ------------//
            std::vector<HalfEdgeEdge*> convexAHalfEdges;

            auto lambda = [&convexAHalfEdges](HalfEdgeEdge* edge)
            {
                convexAHalfEdges.push_back(edge);
            };

            faceA->forEachEdge(lambda);
            facej = 0;


            for (const auto faceB : convexB->GetHalfEdgeFaces())
            {

                //----------------- Get all edges of faceB ------------//
                std::vector<HalfEdgeEdge*> convexBHalfEdges;

                auto lambda = [&convexBHalfEdges](HalfEdgeEdge* edge)
                {
                    convexBHalfEdges.push_back(edge);
                };

                faceB->forEachEdge(lambda);

                bool isCorrectFaces = facei == 0 && facej == 13;

                int edgeAIter = 0;

                for (HalfEdgeEdge* edgeA : convexAHalfEdges)
                {
                    bool drawA = isCorrectFaces && edgeAIter == 1;

                    if (drawA)
                    {
                        //edgeA->DEBUG_drawEdge(transformA, math::colors::magenta, 5.0f);
                    }

                    int edgeBIter = 0;

                    for (HalfEdgeEdge* edgeB : convexBHalfEdges)
                    {
                        bool drawB = isCorrectFaces && edgeBIter == 0;

                        if (drawB)
                        {
                            //edgeB->DEBUG_drawEdge(transformB, math::colors::cyan, 5.0f);
                        }

                        edgeBIter++;

                        //if the given edges creates a minkowski face
                        if (attemptBuildMinkowskiFace(edgeA, edgeB, transformA, transformB))
                        {
                            //get world edge direction
                            math::vec3 edgeADirection = transformA * math::vec4(edgeA->getLocalEdgeDirection(), 0);

                            math::vec3 edgeBDirection = transformB * math::vec4(edgeB->getLocalEdgeDirection(), 0);

                            edgeADirection = math::normalize( edgeADirection );
                            edgeBDirection = math::normalize( edgeBDirection );

                            static float angleThres = math::cos(math::deg2rad(3.0f));
                            float angle = math::abs( math::dot(edgeADirection, edgeBDirection) );

                            if ( angle > angleThres )
                            {
                                continue;
                            }

                            math::vec3 seperatingAxis = math::normalize(math::cross(edgeADirection, edgeBDirection));

                            //get world edge position
                            math::vec3 edgeAtransformedPosition = transformA * math::vec4(edgeA->edgePosition, 1);
                            math::vec3 edgeBtransformedPosition = transformB * math::vec4(edgeB->edgePosition, 1);

                            

                            //check if its pointing in the right direction 
                            if (math::dot(seperatingAxis, edgeAtransformedPosition - positionA) < 0)
                            {
                                seperatingAxis = -seperatingAxis;
                            }


                            //check if given edges create a seperating axis
                            float distance = math::dot(seperatingAxis, edgeBtransformedPosition - edgeAtransformedPosition);
                            //log::debug("distance {} , currentMinimumSeperation {}", distance, currentMinimumSeperation);
                            if (distance > currentMaximumSeperation)
                            {
                                refEdge.ptr = edgeA;
                                incEdge.ptr = edgeB;

                                seperatingAxisFound = seperatingAxis;
                                currentMaximumSeperation = distance;
                            }

                            if (distance > 0)
                            {
                                return true;
                            }

                        }
                        
                    }

                    edgeAIter++;
                }
                facej++;
            }
            facei++;
        }

        maximumSeperation = currentMaximumSeperation;
        return false;
    }

    bool PhysicsStatics::DetectConvexSphereCollision(ConvexCollider* convexA, const math::mat4& transformA, math::vec3 sphereWorldPosition, float sphereRadius,
        float& maximumSeperation)
    {
        //-----------------  check if the seperating axis is the line generated between the centroid of the hull and sphereWorldPosition ------------------//

        math::vec3 worldHullCentroid = transformA * math::vec4(convexA->GetLocalCentroid(), 1);
        math::vec3 centroidSeperatingAxis = math::normalize(worldHullCentroid - sphereWorldPosition);

        math::vec3 seperatingPlanePosition = sphereWorldPosition + centroidSeperatingAxis * sphereRadius;

        math::vec3 worldSupportPoint;
        GetSupportPoint(seperatingPlanePosition, -centroidSeperatingAxis, convexA, transformA, worldSupportPoint);

        float seperation = math::dot(worldSupportPoint - seperatingPlanePosition, centroidSeperatingAxis);

        if (seperation > 0.0f)
        {
            maximumSeperation = seperation;
            return false;
        }

        maximumSeperation = std::numeric_limits<float>::lowest();

        //--------------------------------- check if the seperating axis one of the faces of the convex hull ----------------------------------------------//

        for (auto faceA : convexA->GetHalfEdgeFaces())
        {
            math::vec3 worldFaceCentroid = transformA * math::vec4(faceA->centroid, 1);
            math::vec3 worldFaceNormal = math::normalize(transformA * math::vec4(faceA->normal, 0));

            float seperation = PointDistanceToPlane(worldFaceNormal, worldFaceCentroid, seperatingPlanePosition );

            if (seperation > maximumSeperation)
            {
                maximumSeperation = seperation;
            }

            if (seperation > sphereRadius)
            {
                return false;
            }

        }


        return true;
    }

    std::pair< math::vec3, math::vec3> PhysicsStatics::ConstructAABBFromPhysicsComponentWithTransform
    (ecs::component_handle<physicsComponent> physicsComponentToUse,const math::mat4& transform)
    {
        math::vec3 min, max;

        //auto physicsComponent = physicsComponentToUse.read();

        ////get up
        //math::vec3 invTransUp = math::normalize( math::inverse(transform) * math::vec4(math::vec3(0, 1, 0), 0) );
        //max.y = GetPhysicsComponentSupportPointAtDirection(invTransUp, physicsComponent);
        //
        ////get down
        //math::vec3 invTransDown = math::normalize(math::inverse(transform) * math::vec4(math::vec3(0, -1, 0), 0));
        //min.y = GetPhysicsComponentSupportPointAtDirection(invTransUp, physicsComponent);

        ////get right


        ////get left


        ////get forward


        ////get backward


        return std::make_pair(min,max);
    }

    float PhysicsStatics::GetPhysicsComponentSupportPointAtDirection(math::vec3 direction, physicsComponent& physicsComponentToUse)
    {
        float currentMaximumSupportPoint = std::numeric_limits<float>::lowest();

        //std::vector<math::vec3> vertices;
        ////for each vertex list of each collider
        //for (auto collider : *physicsComponentToUse.colliders)
        //{
        //    auto [first,second] = collider->GetminMaxLocalAABB();
        //    vertices.push_back(first);
        //    vertices.push_back(second);
        //}

        //for (const auto& vert : vertices)
        //{
        //    float dotResult = math::dot(direction, vert);

        //    if (dotResult > currentMaximumSupportPoint)
        //    {
        //        currentMaximumSupportPoint = dotResult;
        //    }
        //}

        return currentMaximumSupportPoint;
    }

    std::pair<math::vec3, math::vec3> PhysicsStatics::ConstructAABBFromVertices(const std::vector<math::vec3>& vertices)
    {
        math::vec3 min, max;

        ////up
        //max.y = GetSupportPoint(vertices, math::vec3(0, 1, 0));
        ////down
        //min.y = GetSupportPoint(vertices, math::vec3(0, -1, 0));

        ////right
        //max.x = GetSupportPoint(vertices, math::vec3(1, 0, 0));
        ////left
        //min.x = GetSupportPoint(vertices, math::vec3(-1, 0, 0));
   

        ////forward
        //max.z = GetSupportPoint(vertices, math::vec3(0, 0, 1));
        ////backward
        //min.z = GetSupportPoint(vertices, math::vec3(0, 0, -1));

        return std::make_pair(min,max);
    }

    std::pair<math::vec3, math::vec3> PhysicsStatics::ConstructAABBFromTransformedVertices(const std::vector<math::vec3>& vertices, const math::mat4& transform)
    {
        math::vec3 min, max;
        math::vec3 worldPos = transform[3];

        math::vec3 outVec;
        //up
        math::vec3 invTransUp = math::normalize(math::inverse(transform) * math::vec4(0, 1, 0, 0));
        GetSupportPoint(vertices, invTransUp, outVec);
        max.y = (transform * math::vec4( outVec,1)).y;

        //down
        math::vec3 invTransDown = math::normalize(math::inverse(transform) * math::vec4(0, -1, 0, 0));
        GetSupportPoint(vertices, invTransDown, outVec);
        min.y = (transform * math::vec4(outVec, 1)).y;

        //right
        math::vec3 invTransRight = math::normalize(math::inverse(transform) * math::vec4(1, 0, 0, 0));
        GetSupportPoint(vertices, invTransRight, outVec);
        max.x = (transform * math::vec4(outVec, 1)).x;

        //left
        math::vec3 invTransLeft = math::normalize(math::inverse(transform) * math::vec4(-1, 0, 0, 0));
        GetSupportPoint(vertices, invTransLeft, outVec);
        min.x = (transform * math::vec4(outVec, 1)).x;

        //forward
        math::vec3 invTransForward = math::normalize(math::inverse(transform) * math::vec4(0, 0, 1, 0));
        GetSupportPoint(vertices, invTransForward, outVec);
        max.z = (transform * math::vec4(outVec, 1)).z;

        //backward
        math::vec3 invTransBackward = math::normalize(math::inverse(transform) * math::vec4(0, 0, -1, 0));
        GetSupportPoint(vertices, invTransBackward, outVec);
        min.z = (transform * math::vec4(outVec, 1)).z;

      

        return std::make_pair(min, max);
    }

    std::pair<math::vec3, math::vec3> PhysicsStatics::CombineAABB(const std::pair<math::vec3, math::vec3>& first, const std::pair<math::vec3, math::vec3>& second)
    {
        auto& firstLow = first.first;
        auto& firstHigh = first.second;
        auto& secondLow = second.first;
        auto& secondHigh = second.second;
        math::vec3 lowBounds = secondLow;
        math::vec3 highBounds = secondHigh;
        if (firstLow.x < secondLow.x)   lowBounds.x    = firstLow.x;
        if (firstLow.y < secondLow.y)   lowBounds.y    = firstLow.y;
        if (firstLow.z < secondLow.z)   lowBounds.z    = firstLow.z;
        if (firstHigh.x > secondHigh.x) highBounds.x   = firstHigh.x;
        if (firstHigh.y > secondHigh.y) highBounds.y   = firstHigh.y;
        if (firstHigh.z > secondHigh.z) highBounds.z   = firstHigh.z;

        return std::make_pair(lowBounds, highBounds);
    }

    float PhysicsStatics::FindClosestPointToLineInterpolant(const math::vec3& startPoint, const math::vec3& lineDirection, const math::vec3& pointPosition)
    {
        return ( math::dot(lineDirection,pointPosition) - math::dot(lineDirection,startPoint) ) / math::dot(lineDirection,lineDirection);
    }

    math::vec3 PhysicsStatics::FindClosestPointToLineSegment(const math::vec3& start, const math::vec3& end, const math::vec3& pointPosition)
    {
        float interpolant = FindClosestPointToLineInterpolant(start, end - start, pointPosition);
        interpolant = math::clamp(interpolant, 0.0f, 1.0f);

        return start + (end - start) * interpolant;
    }

    std::shared_ptr<ConvexCollider> PhysicsStatics::GenerateConvexHull(const std::vector<math::vec3>& vertices, int maxDraw, int DEBUG_at, math::mat4 DEBUG_transform)
    {
        //[1] Calculate scaled epsilon
        static float initialEpsilon = math::sqrt(math::epsilon<float>());

        math::vec3 maxInDimension(std::numeric_limits<float>::lowest());

        for (const auto& vert : vertices)
        {
            if (math::abs(vert.x) > maxInDimension.x)
            {
                maxInDimension.x = math::abs(vert.x);
            }

            if (math::abs(vert.y) > maxInDimension.y)
            {
                maxInDimension.y = math::abs(vert.y);
            }

            if (math::abs(vert.z) > maxInDimension.z)
            {
                maxInDimension.z = math::abs(vert.z);
            }
        }
        
        //epsilon must take into account span of vertices
        float dimensionSum = 3 * (math::abs(maxInDimension.x) + math::abs(maxInDimension.y) + math::abs(maxInDimension.z));
        const float scaledEpsilon = dimensionSum * initialEpsilon;

        std::vector<HalfEdgeFace*> faces;
        faces.reserve(4);

        //[3] Measure extents of mesh
        std::array<math::vec3, 6> supportVertices;

        //(1.1) Get support points in -x and x
        GetSupportPoint(vertices, math::vec3(1, 0, 0), supportVertices.at(0));
        GetSupportPoint(vertices, math::vec3(-1, 0, 0), supportVertices.at(1));

        //(1.2) Get support points in -y and y
        GetSupportPoint(vertices, math::vec3(0, 1, 0), supportVertices.at(2));
        GetSupportPoint(vertices, math::vec3(0, -1, 0), supportVertices.at(3));

        //(1.3) Get support points in -z and z
        GetSupportPoint(vertices, math::vec3(0, 0, 1), supportVertices.at(4));
        GetSupportPoint(vertices, math::vec3(0, 0, -1), supportVertices.at(5));

        float x = math::abs(supportVertices.at(0).x - supportVertices.at(1).x);
        float y = math::abs(supportVertices.at(2).y - supportVertices.at(3).y);
        float z = math::abs(supportVertices.at(4).z - supportVertices.at(5).z);

        math::vec3 start ( supportVertices.at(1).x, supportVertices.at(3).y, supportVertices.at(5).z);
        start = DEBUG_transform * math::vec4(start, 1);

       /* debug::drawLine(start, start + math::vec3(x, 0, 0), math::colors::red, 20.0f, FLT_MAX, true);
        debug::drawLine(start, start + math::vec3(0, y, 0), math::colors::blue, 20.0f, FLT_MAX, true);
        debug::drawLine(start, start + math::vec3(0, 0, z), math::colors::green, 20.0f, FLT_MAX, true);*/

        //merge volume threshold is equal a certain small fraction of the volume of the bounding box 
        float mergeVolumeThreshold = x * y * z * (1.0f/3000.0f);

        //[2] Build Initial Hull
        if (!qHBuildInitialHull(vertices, supportVertices, faces,DEBUG_transform))
        {
            return nullptr;
        }

        int currentDraw = 0;
        

        //[3] populate list with faces of initial hull
        std::list<ColliderFaceToVert> facesWithOutsideVerts;
        partitionVerticesToList(vertices, faces, facesWithOutsideVerts,false);

        //[4] loop through faces until there are no faces with unmerged vertices
        if (!facesWithOutsideVerts.empty())
        {
            PointerEncapsulator< ColliderFaceToVert> currentFaceToVert;
            currentFaceToVert.ptr = &facesWithOutsideVerts.front();
   
            while (foundFaceWithOutsideVert(facesWithOutsideVerts, currentFaceToVert))
            {
                //log::debug("Iter {0}", currentDraw);

               /* if (currentDraw == 29)
                {
                    DebugBreak();
                }*/

                //if (currentDraw >= maxDraw) { break; }
                currentDraw++;

                bool atDebug = currentDraw == DEBUG_at + 1;

                //find furhtest vertex of last face
                auto [furthestVert, distanceFromFace] = currentFaceToVert.ptr->GetFurthestOutsideVert();
                math::vec3 worldPos = DEBUG_transform * math::vec4(furthestVert, 1);

                //debug::drawLine(worldPos, worldPos + math::vec3(0, 0.1, 0), math::colors::magenta, 5.0f, FLT_MAX, false);

                //check if we should merge this vertex
                if (distanceFromFace > scaledEpsilon)
                {
                    bool sucess = mergeVertexToHull(furthestVert, facesWithOutsideVerts,
                        scaledEpsilon, mergeVolumeThreshold, DEBUG_transform, atDebug);

                    if (!sucess)
                    {
                        currentFaceToVert.ptr->outsideVerts.clear();
                    }
                }
                else
                {
                    //this face has no mergable vertices
                    currentFaceToVert.ptr->outsideVerts.clear();
                }

                if (currentDraw == DEBUG_at)
                {
                    PointerEncapsulator< ColliderFaceToVert> debugFaceToVert;
                    foundFaceWithOutsideVert(facesWithOutsideVerts, debugFaceToVert);
                    auto [furthestVert, distanceFromFace] = debugFaceToVert.ptr->GetFurthestOutsideVert();
                    math::vec3 worldPos = DEBUG_transform * math::vec4(furthestVert, 1);
                    debug::drawLine(worldPos, worldPos + math::vec3(0, 0.05f, 0), math::colors::magenta, 5.0f, FLT_MAX, true);
                }

                if (atDebug)
                {
                    break;
                }
            }

        }


        auto convexCollider = std::make_shared<ConvexCollider>();
        auto& halfEdgesVector = convexCollider->GetHalfEdgeFaces();

        for (auto& faceToVert : facesWithOutsideVerts)
        {
            
            halfEdgesVector.push_back(faceToVert.face);
        }
        

        //populate list of vertices in collider list

        convexCollider->PopulateVertexListWithHalfEdges();

        return convexCollider;
    }

    void PhysicsStatics::CalculateNewellPlane(const std::vector<math::vec3>& v, math::vec3& outPlaneNormal, float& distToCentroid)
    {
        math::vec3 centroid{0,0,0};
        outPlaneNormal = math::vec3();

        for (int i = v.size() - 1, j = 0; j < v.size(); i = j, j++)
        {
            outPlaneNormal.x += (v[i].y - v[j].y) * (v[i].z + v[j].z); // projection on yz
            outPlaneNormal.y += (v[i].z - v[j].z) * (v[i].x + v[j].x); // projection on xz
            outPlaneNormal.z += (v[i].x - v[j].x) * (v[i].y + v[j].y); // projection on xy
            centroid += v[j];
        }
           
        // Normalize normal and fill in the plane equation fields
        outPlaneNormal = math::normalize(outPlaneNormal);
        distToCentroid = math::dot(centroid, outPlaneNormal) / v.size(); // “centroid / n” is the true centroid point
    }

    bool PhysicsStatics::attemptBuildMinkowskiFace(HalfEdgeEdge* edgeA, HalfEdgeEdge* edgeB, const math::mat4& transformA, const math::mat4& transformB)
    {
        //TODO the commmented parts are technically more robust and should work but somehow dont, figure out why.
        //11 , 745
        const math::vec3 transformedA1 = transformA *
            math::vec4(edgeA->getLocalNormal(), 0);

        const math::vec3 transformedA2 = transformA *
            math::vec4(edgeA->pairingEdge->getLocalNormal(), 0);

        const math::vec3 transformedEdgeDirectionA = transformA * math::vec4(edgeA->getLocalEdgeDirection(), 0);
        //

        const math::vec3 transformedB1 = transformB *
            math::vec4(edgeB->getLocalNormal(), 0);

        const math::vec3 transformedB2 = transformB *
            math::vec4(edgeB->pairingEdge->getLocalNormal(), 0);

        const math::vec3 transformedEdgeDirectionB = transformB * math::vec4(edgeB->getLocalEdgeDirection(), 0);
        //
        //if (edgeA->identifier == 49 && edgeB->identifier == 484)
        //{
        //    math::vec3 positionA = transformA * math::vec4(edgeA->edgePosition, 1);
        //    debug::drawLine(positionA, positionA + transformedA1, math::colors::green, 5.0f, FLT_MAX, true);
        //    debug::drawLine(positionA, positionA + transformedA2, math::colors::green, 5.0f, FLT_MAX, true);

        //    debug::drawLine(positionA, positionA - transformedB1, math::colors::red, 5.0f, FLT_MAX, true);
        //    debug::drawLine(positionA, positionA - transformedB2, math::colors::red, 5.0f, FLT_MAX, true);
        //    auto x = 1;
        //}

        return isMinkowskiFace(transformedA1, transformedA2, -transformedB1, -transformedB2
            , transformedEdgeDirectionA, transformedEdgeDirectionB);
    }

    bool PhysicsStatics::isMinkowskiFace(const math::vec3& transformedA1, const math::vec3& transformedA2, const math::vec3& transformedB1,
        const math::vec3& transformedB2, const math::vec3& planeANormal, const math::vec3& planeBNormal)
    {
        //------------------------ Check if normals created by arcA seperate normals of B --------------------------------------//
        //CBA
        float planeADotB1 = math::dot(planeANormal, transformedB1);
        //DBA
        float planeADotB2 = math::dot(planeANormal, transformedB2);

        float dotMultiplyResultA =
            planeADotB1 * planeADotB2;

        //log::debug("dotMultiplyResultA {}", dotMultiplyResultA);

        if (dotMultiplyResultA > 0.0f )
        {
            return false;
        }

        //------------------------ Check if normals created by arcB seperate normals of A --------------------------------------//

        //ADC
        float planeBDotA1 = math::dot(planeBNormal, transformedA1);
        //BDC
        float planeBDotA2 = math::dot(planeBNormal, transformedA2);

        float  dotMultiplyResultB = planeBDotA1 * planeBDotA2;

        //log::debug("dotMultiplyResultB {}", dotMultiplyResultB);

        if (dotMultiplyResultB > 0.0f )
        {
            return false;
        }

        //------------------------ Check if arcA and arcB are in the same hemisphere --------------------------------------------//

        math::vec3 abNormal = math::cross(transformedA2, transformedB2);

        float planeABDotA1 = math::dot(abNormal, transformedA1);
        float planeABDotB1 = math::dot(abNormal, transformedB1);

        float dotMultiplyResultAB = planeABDotA1 * planeABDotB1;

        if (dotMultiplyResultAB <= 0.0f )
        {
            return false;
        }

        return true;
    }

    bool PhysicsStatics::qHBuildInitialHull(const std::vector<math::vec3>& vertices,
        std::array<math::vec3,6>& supportVertices, std::vector<HalfEdgeFace*>& faces, math::mat4 DEBUG_transform )
    {
        //Summary:
        //[1] Find the 2 most distant vertices in 'support Vertices'
        //[2] Find the vertex most distant from the line created by the 2 most distance vertices
        //[3] Create first collider face using that line and the vertex most distant to it
        //[4] Find the most distant vertex from the plane where the first collider face lies
        //[5] Create a set of faces connecting the first collider face to the most distant vertex,
        //this create a tetrahedron shaped collider

        //[1] Find the 2 most distant vertices in 'vertices'

        float mostDistant = std::numeric_limits<float>::lowest();

        math::vec3& firstDistant = supportVertices.at(0);
        math::vec3& secondDistant = supportVertices.at(1);

        //(1.1) Iterate through support vertices to find the combination of vertices that represent the 2 most distant points
        for (int i = 0; i < supportVertices.size(); i++)
        {
            math::vec3& first = supportVertices.at(i);

            for (int j = i+1; j < supportVertices.size()-1; j++)
            {
                math::vec3& second = supportVertices.at(j);

                float currentDistance2 = math::distance2(first, second);

                if (currentDistance2 > mostDistant)
                {
                    mostDistant = currentDistance2;
                    firstDistant = first;
                    secondDistant = second;
                }
            }
        }

        //DEBUG draw most distant
      /*  {
            math::vec3 transformFirst = DEBUG_transform * math::vec4(firstDistant, 1);
            math::vec3 transformSecond = DEBUG_transform * math::vec4(secondDistant, 1);

            debug::drawLine(transformFirst, transformSecond,math::colors::red,5.0f,FLT_MAX,true);
        }*/

        //[2] Find the vertex most distant from the line created by the 2 most distance vertices

        math::vec3 firstToSecond = math::normalize(secondDistant - firstDistant);

        const math::vec3* thirdDistant = nullptr;
        mostDistant = std::numeric_limits<float>::lowest();

        //(2.1) Iterate through 'vertices' to find the vertex most distant from line 
        for (auto& vertex : vertices)
        {
            //(2.1.1) Check if vertex and firstDistant create a line segment parralel to the line segment created by 'firstToSecond'
            float dotResult = math::dot(math::normalize(vertex - firstDistant), firstToSecond);

            if (math::close_enough(dotResult, 1.0f))
            {
                continue;
            }

            //(2.1.2) Find closest point between vertex and line segment
            math::vec3 closestPoint = PhysicsStatics::FindClosestPointToLineSegment(firstDistant, secondDistant, vertex);

            float currentDistance = math::distance2(closestPoint, vertex);

            //(2.1.3) Store vertex if it is further than the current known most distant
            if (currentDistance > mostDistant)
            {
                mostDistant = currentDistance;
                thirdDistant = &vertex;
            }
        }

        if (!thirdDistant)
        {
            return false;
        }

        //DEBUG draw most distant
      /*  {
            math::vec3 transformFirst = DEBUG_transform * math::vec4(*thirdDistant, 1);

            debug::drawLine(transformFirst, transformFirst + math::vec3(0,0.1f,0), math::colors::blue, 5.0f, FLT_MAX, true);
        }*/

        //[3] Create first collider face using that line and the vertex most distant to it

        //(3.1) Initialize Half Edges
        HalfEdgeEdge* firstEdge = new HalfEdgeEdge(firstDistant);
        HalfEdgeEdge* secondEdge = new HalfEdgeEdge(secondDistant);
        HalfEdgeEdge* thirdEdge = new HalfEdgeEdge(*thirdDistant);

        //(3.2) Connect them to each other
        firstEdge->setNextAndPrevEdge(thirdEdge, secondEdge);
        secondEdge->setNextAndPrevEdge(firstEdge, thirdEdge);
        thirdEdge->setNextAndPrevEdge(secondEdge, firstEdge);

        //(3.3) Initialize Half Edge Faces
        HalfEdgeFace* initialFace = new HalfEdgeFace(firstEdge,
            math::normalize(math::cross(secondDistant - firstDistant, *thirdDistant - secondDistant)));

        //(3.4) Add to collider
        faces.push_back(initialFace);

        //[4] Find the most distant vertex from the plane where the first collider face lies

        //(4.1) Iterate through vertices to the a point to plane check
        mostDistant = std::numeric_limits<float>::lowest();
        const math::vec3* firstEyePoint = nullptr;

        math::vec3 planePosition = initialFace->centroid;
        math::vec3 planeNormal = initialFace->normal;

        for (auto& vertex : vertices)
        {
            float currentDistance =
                math::abs(PhysicsStatics::PointDistanceToPlane(planeNormal, planePosition, vertex));

            if (currentDistance > mostDistant)
            {
                mostDistant = currentDistance;
                firstEyePoint = &vertex;
            }
        }

        if (math::close_enough(mostDistant,0.0f))
        {
            return false;
        }

        math::color drawC = math::colors::blue;
       

        //[5] invert face if distant vertex is in front of face

        float eyePointDistance =
            PhysicsStatics::PointDistanceToPlane(planeNormal, planePosition, *firstEyePoint);
        bool needInverse = eyePointDistance > 0.0f;
        
        if (needInverse)
        {
            initialFace->inverse();
            drawC = math::colors::red;
        }

        //DEBUG draw most distant
        /*{
            math::vec3 transformFirst = DEBUG_transform * math::vec4(*firstEyePoint, 1);

            debug::drawLine(transformFirst, transformFirst + math::vec3(0, 0.1f, 0), drawC, 5.0f, FLT_MAX, true);

        }*/

        //[5] Create a set of faces connecting the first collider face to the most distant vertex

        //(5.1) Reverse collect the edges of the initialFace

        std::vector<HalfEdgeEdge*> reverseHalfEdgeList;

        auto collectEdges = [&reverseHalfEdgeList](HalfEdgeEdge* current)
        {
            reverseHalfEdgeList.push_back(current);
        };

        initialFace->forEachEdgeReverse(collectEdges);

        //(5.2) For each edge create a new face that connects the initial face with the eyePoint
        math::vec3 eyePoint = *firstEyePoint;

        createHalfEdgeFaceFromEyePoint(eyePoint, reverseHalfEdgeList,faces);

        return true;
    }

    void PhysicsStatics::createHalfEdgeFaceFromEyePoint(const math::vec3 eyePoint, const std::vector<HalfEdgeEdge*>& reversedEdges, std::vector<HalfEdgeFace*>& createdFaces)
    {
        HalfEdgeEdge* pairingToConnectTo = nullptr;
        HalfEdgeEdge* initialPairing = nullptr;

        for (auto edge : reversedEdges)
        {
            //initialize pairing its position is on next
            HalfEdgeEdge* pairing = new HalfEdgeEdge(edge->nextEdge->edgePosition);

            //initialize next pairing its position is on current
            HalfEdgeEdge* nextPairing = new HalfEdgeEdge(edge->edgePosition);

            //initialize prev pairing its position is on the eye point
            HalfEdgeEdge* prevPairing = new HalfEdgeEdge(eyePoint);

            //connect to each other
            pairing->setNextAndPrevEdge(prevPairing, nextPairing);
            nextPairing->setNextAndPrevEdge(pairing, prevPairing);
            prevPairing->setNextAndPrevEdge(nextPairing, pairing);

            //
            pairing->setPairingEdge(edge);

            //initialize new face
            math::vec3 faceNormal = math::normalize(math::cross(nextPairing->edgePosition - pairing->edgePosition,
                prevPairing->edgePosition - pairing->edgePosition));

            HalfEdgeFace* face = new HalfEdgeFace(pairing, faceNormal);

            //push new face into list 
            createdFaces.push_back(face);

            //connect to 
            if (pairingToConnectTo)
            {
                prevPairing->setPairingEdge(pairingToConnectTo);
            }
            else
            {
                initialPairing = prevPairing;
            }

            pairingToConnectTo = nextPairing;
        }

        initialPairing->setPairingEdge(pairingToConnectTo);
    }

    bool PhysicsStatics::foundFaceWithOutsideVert(std::list<ColliderFaceToVert>& facesWithOutsideVerts, PointerEncapsulator< ColliderFaceToVert>& outChosenFace)
    {
        int i = 0;
        for (auto& faceWithOutsideVert : facesWithOutsideVerts)
        {
            if (!faceWithOutsideVert.outsideVerts.empty())
            {
                outChosenFace.ptr = &faceWithOutsideVert;
                return true;
            }
            i++;
        }

        return false;
    }

    void PhysicsStatics::partitionVerticesToList(const std::vector<math::vec3> vertices, const std::vector<HalfEdgeFace*>& faces,
        std::list<ColliderFaceToVert>& outFacesWithOutsideVerts, bool shouldFindFirstMatch )
    {

        auto findFirstMatch = [&outFacesWithOutsideVerts](const math::vec3& vertex)
        {
            for (ColliderFaceToVert& faceToVert : outFacesWithOutsideVerts)
            {
                if (IsPointAbovePlane(faceToVert.face->normal, faceToVert.face->centroid, vertex))
                {
                    faceToVert.outsideVerts.push_back(vertex);
                    return;
                }
            }
        };

        auto findBestMatch = [&outFacesWithOutsideVerts]( const math::vec3& vertex)
        {
            ColliderFaceToVert* bestFaceToVert = nullptr;
            float bestMatchDistance = std::numeric_limits<float>::lowest();
           

            for (ColliderFaceToVert& faceToVert : outFacesWithOutsideVerts)
            {
                auto face = faceToVert.face;
                float currentDistance = PointDistanceToPlane(face->normal, face->centroid, vertex);

                if (currentDistance > bestMatchDistance)
                {
                    bestMatchDistance = currentDistance;
                    bestFaceToVert = &faceToVert;

                }
            }

            if (bestFaceToVert && bestMatchDistance > 0.0f)
            {
                bestFaceToVert->outsideVerts.push_back(vertex);
            }
        };

        legion::core::delegate<void(const math::vec3& vertex)> partioningStrategyFunc = nullptr;

        if (shouldFindFirstMatch)
        {
            partioningStrategyFunc = findFirstMatch;
        }
        else
        {
            partioningStrategyFunc = findBestMatch;
        }

        for (HalfEdgeFace* face : faces)
        {
            outFacesWithOutsideVerts.emplace_back(face);
        }

        //for each vertex in vertices
        for (const math::vec3& vertex : vertices)
        {
            partioningStrategyFunc(vertex);
        }
    }

    void PhysicsStatics::findHorizonEdgesFromFaces(const math::vec3& eyePoint, std::vector<HalfEdgeFace*>& faces,
        std::vector<HalfEdgeEdge*>& outHorizonEdges, float scalingEpsilon, math::mat4 DEBUG_transform, bool atDebug)
    {
        //[1] Find first horizon edge
        HalfEdgeEdge* initialHorizon = nullptr;

        auto findFirstHorizon = [&eyePoint,&initialHorizon,scalingEpsilon](HalfEdgeEdge*edge)
        {
            if (initialHorizon) { return; }
           
            if (edge->isEdgeHorizonFromVertex(eyePoint, scalingEpsilon))
            {
                initialHorizon = edge;
            }
        };

        for (auto face : faces)
        {
            if (initialHorizon) { break; }

            face->forEachEdge(findFirstHorizon);
            
        }

        assert(initialHorizon);

        //[2] Loop through the collider to collect the other horizon edges
       
        HalfEdgeEdge* currentEdge = initialHorizon;
        //outHorizonEdges.push_back(currentEdge);
        if (atDebug)
        {
            initialHorizon->DEBUG_drawEdge(DEBUG_transform, math::colors::orange, FLT_MAX, 5.0f);
            //return;
        }
        
        int iter = 0;
        do
        {

            currentEdge = currentEdge->nextEdge;

            if (atDebug) { log::debug("-> Iteration"); currentEdge->DEBUG_drawEdge(DEBUG_transform, math::colors::cyan, FLT_MAX, 5.0f); }

            if (!currentEdge->isEdgeHorizonFromVertex(eyePoint, scalingEpsilon))
            {
                if (atDebug) { log::debug("Next edge was not a horizon, finding horizon now"); }
                int safetyCount = 0;
                auto stuckEdge = currentEdge;
                do
                {
                    if (atDebug) { log::debug("getting pairingEdge->nextEdge"); }
                    auto pairing = currentEdge->pairingEdge;
                    assert(pairing != currentEdge);
                    currentEdge = currentEdge->pairingEdge->nextEdge;
               
                    if (iter == 2 && atDebug)
                    {
                        pairing->DEBUG_directionDrawEdge(DEBUG_transform, math::colors::blue, FLT_MAX, 40.0f);
                        currentEdge->DEBUG_directionDrawEdge(DEBUG_transform, math::colors::blue, FLT_MAX, 40.0f);
                        //currentEdge->pairingEdge->DEBUG_directionDrawEdge(DEBUG_transform, math::colors::blue, FLT_MAX, 40.0f);
                        //DebugBreak();
                        return;
                    }
                    assert(currentEdge);
                   /* safetyCount++;

                    if (safetyCount > 10)
                    {
                        stuckEdge->DEBUG_drawEdge(DEBUG_transform, math::colors::red, FLT_MAX, 5.0f);
                        return;
                    }*/
                    
                }
                while (!currentEdge->isEdgeHorizonFromVertex(eyePoint, scalingEpsilon));

            }

            //if (atDebug) { log::debug("getting nextEdge"); }

            outHorizonEdges.push_back(currentEdge);
            iter++;

        } while (currentEdge != initialHorizon);


        for (auto edge : outHorizonEdges)
        {
            //edge->DEBUG_drawEdge(DEBUG_transform, math::colors::red, FLT_MAX, 5.0f);
        }

    }

    bool PhysicsStatics::mergeVertexToHull(const math::vec3& eyePoint,std::list<ColliderFaceToVert>& facesWithOutsideVerts,
        float scalingEpsilon, float hullMinimumVolume, math::mat4 DEBUG_transform, bool atDebug)
    {
        //if (atDebug) { return true; }
        auto testPairings = [](HalfEdgeEdge* edge)
        {
            math::vec3 x = edge->pairingEdge->edgePosition;
        };

        std::vector<math::vec3> unmergedVertices;
        std::vector<HalfEdgeFace*> facesToBeRemoved;

        //TODO we are doing the point to distance thing twice here, would be better if we reused the results somehow
        float totalVolume = 0.0f;
        std::vector<HalfEdgeFace*> visibleFaces;
        //[1] Calculate volume created if this particular vertex was merged
        for (auto listIter = facesWithOutsideVerts.begin(); listIter != facesWithOutsideVerts.end(); listIter++)
        {
            HalfEdgeFace* face = listIter->face;

            const math::vec3& planeCentroid = face->centroid;
            const math::vec3& planeNormal = face->normal;

            float distanceToPlane = PointDistanceToPlane(planeNormal, planeCentroid, eyePoint);

            if (distanceToPlane > scalingEpsilon)
            {
                float faceArea = listIter->face->CalculateFaceArea();

                totalVolume += faceArea * distanceToPlane * (1.0f/3.0f);
                visibleFaces.push_back(face);
            }
        }

        //[2] Only continue if volume is above a given threshold
        if (totalVolume < hullMinimumVolume)
        {
            static bool notFirstCancelledHull = true;

            if (notFirstCancelledHull)
            {
                notFirstCancelledHull = false;
                //log::debug("First Cancelled HUll");

                /*for (auto face : visibleFaces)
                {
                    face->DEBUG_DrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
                }

                math::vec3 worldPos = DEBUG_transform * math::vec4(eyePoint, 1);
                debug::drawLine(worldPos, worldPos + math::vec3(0, 0.1, 0), math::colors::magenta, 5.0f, FLT_MAX, true);*/

            }

            return false;
        }


        //identify faces that can see vertex and remove them from list
        for (auto listIter = facesWithOutsideVerts.begin(); listIter != facesWithOutsideVerts.end();)
        {
            HalfEdgeFace* face = listIter->face;

            const math::vec3& planeCentroid = face->centroid;
            const math::vec3& planeNormal = face->normal;

            float distanceToPlane = PointDistanceToPlane(planeNormal, planeCentroid, eyePoint);

            if (distanceToPlane > scalingEpsilon)
            {
                //face can see vertex, we must remove it from list
                facesToBeRemoved.push_back(face);
                listIter->populateVectorWithVerts(unmergedVertices);
                listIter = facesWithOutsideVerts.erase(listIter);
            }
            else
            {
                listIter++;
            }
        }

        if (atDebug)
        {
            for (auto face : facesToBeRemoved)
            {
                //face->DEBUG_DirectionDrawFace(DEBUG_transform, math::colors::grey, FLT_MAX);
            }
            log::debug(" facesToBeRemoved {0} ", facesToBeRemoved.size());
            //return true;
        }
        
        //identify horizon edges and put them into list
        std::vector<HalfEdgeEdge*> horizonEdges;
        findHorizonEdgesFromFaces(eyePoint, facesToBeRemoved, horizonEdges,scalingEpsilon, DEBUG_transform,false);

        if (atDebug)
        {
            //return;
            for (auto edge : horizonEdges)
            {
                //edge->DEBUG_drawEdge(DEBUG_transform, math::colors::red, FLT_MAX);
            }
            log::debug("at debug horizon edges {0} ", horizonEdges.size());
            //return true;
        }

        //reverse iterate the list to find their pairings, add them to new list
        {
            std::vector<HalfEdgeEdge*> tempEdges = std::move(horizonEdges);
            horizonEdges.clear();

            //what we actually want is its pairings
            for (auto horizonEdge : tempEdges)
            {
                horizonEdges.push_back(horizonEdge->pairingEdge);
            }
        }

        if (atDebug)
        {
           /* for (auto face : facesToBeRemoved)
            {
                face->DEBUG_DirectionDrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
            }
            log::debug(" facesToBeRemoved {0} ", facesToBeRemoved.size());*/
            //return true;
        }


        std::vector<HalfEdgeFace*> newFaces;
        newFaces.reserve(horizonEdges.size());
        createHalfEdgeFaceFromEyePoint(eyePoint, horizonEdges, newFaces);

        if (atDebug)
        {
            //433type 
            /*for (size_t i = 0; i < newFaces.size(); i++)
            {
                newFaces.at(i)->DEBUG_DrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
            }*/
            //newFaces.at(0)->DEBUG_DrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
            //newFaces.at(1)->DEBUG_DrawFace(DEBUG_transform, math::colors::blue, FLT_MAX);
            //newFaces.at(2)->DEBUG_DrawFace(DEBUG_transform, math::colors::red, FLT_MAX);

            log::debug(" newFaces {0} ", newFaces.size());
            //return true;
        }
        
        for (int i = 0; i < horizonEdges.size(); i++)
        {
            //horizon is now 
            HalfEdgeFace* establishedFace = horizonEdges.at(i)->face;
            ColliderFaceToVert& faceToVertEstablished = *establishedFace->faceToVert;
            HalfEdgeFace* newFace = newFaces.at(i);

            if (atDebug && i ==1)
            {
                horizonEdges.at(i)->DEBUG_directionDrawEdge(DEBUG_transform, math::colors::blue, FLT_MAX, 4.0f);
                establishedFace->DEBUG_DrawFace(DEBUG_transform, math::colors::green, FLT_MAX);
                newFace->DEBUG_DrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
            }

            if (isNewellFacesCoplanar(establishedFace, newFace, horizonEdges.at(i)->pairingEdge,scalingEpsilon,DEBUG_transform,atDebug))
            {
                //assert(establishedFace != horizonEdges.at(i)->face);
                horizonEdges.at(i)->suicidalMergeWithPairing(DEBUG_transform);

                faceToVertEstablished.populateVectorWithVerts(unmergedVertices);
                newFaces.at(i) = nullptr;
                horizonEdges.at(i) = nullptr;

            }
        }

        //for each edge in horizon edge
            //if edge is available
        for (size_t i = 0; i < horizonEdges.size(); i++)
        {
            size_t next = (i+1) % horizonEdges.size();

            auto nextEdge = horizonEdges.at(next);

            if (newFaces.at(i) && newFaces.at(next))
            {
                HalfEdgeEdge* connectingEdgeFromSecond = nextEdge->pairingEdge->prevEdge;

                /*if (atDebug)
                {
                    connectingEdgeFromSecond->DEBUG_directionDrawEdge(DEBUG_transform, math::colors::magenta, FLT_MAX, 4.0f);
                    nextEdge->DEBUG_directionDrawEdge(DEBUG_transform, math::colors::blue, FLT_MAX, 4.0f);
                    
                    log::debug("-> merged");
                }*/
              
                if (isNewellFacesCoplanar(newFaces.at(i), newFaces.at(next), connectingEdgeFromSecond, scalingEpsilon, DEBUG_transform, atDebug))
                {
                    connectingEdgeFromSecond->suicidalMergeWithPairing(DEBUG_transform);
                    newFaces.at(i) = nullptr;
                }
            }
        }


        {
            std::vector<HalfEdgeFace*> tempNewFaces = std::move(newFaces);

            for (auto face : tempNewFaces)
            {
                if (face)
                {
                    newFaces.push_back(face);
                }
            }
        }

        if (atDebug)
        {
            ////433
            for (auto face : newFaces)
            {
                //face->DEBUG_DrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
            }

           /* newFaces.at(0)->DEBUG_DrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
            newFaces.at(0)->forEachEdge(testPairings);*/
            //newFaces.at(1)->DEBUG_DirectionDrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
            //newFaces.at(1)->forEachEdge(testPairings); 
            //newFaces.at(1)->DEBUG_DrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
           /* newFaces.at(0)->DEBUG_DrawFace(DEBUG_transform, math::colors::red, FLT_MAX);
            newFaces.at(1)->DEBUG_DrawFace(DEBUG_transform, math::colors::blue, FLT_MAX);*/
            log::debug(" final newFaces {0} ", newFaces.size());
            //return true;
        }

        partitionVerticesToList(unmergedVertices, newFaces, facesWithOutsideVerts,false);

        for (auto face : facesToBeRemoved)
        {
            delete face;
        }

        return true;
    }

    bool PhysicsStatics::isFacesCoplanar(HalfEdgeFace* first, HalfEdgeFace* second)
    {
        static float angle = math::cos(math::radians(10.0f));

        float dotResult = math::dot(math::normalize(first->normal), math::normalize(second->normal));
        return dotResult > angle;
    }

    bool PhysicsStatics::isNewellFacesCoplanar(HalfEdgeFace* first, HalfEdgeFace* second,
        HalfEdgeEdge* connectingEdge, float scalingEpsilon, math::mat4 DEBUG_transform, bool atDebug)
    {
        //each face will at least have 3 vertices
        std::vector<math::vec3> NewellPolygon;
        NewellPolygon.reserve(6);

        auto firstOriginal = first->startEdge;
        auto secondOriginal = second->startEdge;

        first->startEdge = connectingEdge->pairingEdge->nextEdge;
        second->startEdge = connectingEdge->nextEdge;

        auto collectVerticesOfFace = [&NewellPolygon](HalfEdgeEdge* edge)
        {
            NewellPolygon.emplace_back(edge->edgePosition);
        };

        first->forEachEdge(collectVerticesOfFace);
        NewellPolygon.pop_back();

        second->forEachEdge(collectVerticesOfFace);
        NewellPolygon.pop_back();

       /* first->startEdge = firstOriginal;
        second->startEdge = secondOriginal;*/


        math::vec3 planeNormal; float distToCentroid;
        CalculateNewellPlane(NewellPolygon, planeNormal, distToCentroid);

        if (atDebug)
        {
            math::vec3 pos;

           /* pos = DEBUG_transform * math::vec4(NewellPolygon.at(0), 1);
            debug::drawLine(pos, pos + math::vec3(0, 0.1f, 0), math::colors::red, 4.0f, FLT_MAX, true);

            pos = DEBUG_transform * math::vec4(NewellPolygon.at(1), 1);
            debug::drawLine(pos, pos + math::vec3(0, 0.1f, 0), math::colors::green, 4.0f, FLT_MAX, true);

            pos = DEBUG_transform * math::vec4(NewellPolygon.at(2), 1);
            debug::drawLine(pos, pos + math::vec3(0, 0.1f, 0), math::colors::blue, 4.0f, FLT_MAX, true);

            pos = DEBUG_transform * math::vec4(NewellPolygon.at(3), 1);
            debug::drawLine(pos, pos + math::vec3(0, 0.1f, 0), math::colors::magenta, 4.0f, FLT_MAX, true);

            pos = DEBUG_transform * math::vec4(NewellPolygon.at(4), 1);
            debug::drawLine(pos, pos + math::vec3(0, 0.1f, 0), math::colors::cyan, 4.0f, FLT_MAX, true);*/

        }


        for (int i = 0; i < NewellPolygon.size(); i++) {
            float dist = math::dot(planeNormal, NewellPolygon[i]) - distToCentroid;

            if (math::abs(dist) > scalingEpsilon)
            {
                return false;
            }
        }

        return true;
    }

};

