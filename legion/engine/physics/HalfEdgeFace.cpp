#include <physics/halfedgeface.hpp>
#include <rendering/debugrendering.hpp>
#include <physics/physics_statics.hpp>

namespace legion::physics
{
    HalfEdgeFace::HalfEdgeFace(HalfEdgeEdge* newStartEdge, math::vec3 newNormal) : startEdge{ newStartEdge }, normal{ newNormal }
    {
        initializeFace();
    }

    void HalfEdgeFace::initializeFace()
    {
        /*log::debug("HalfEdgeFace::HalfEdgeFace");*/
        static int faceCount = 0;

        DEBUG_color = math::color(math::linearRand(0.25f, 0.7f), math::linearRand(0.25f, 0.7f), math::linearRand(0.25f, 0.7f));
        math::vec3 faceCenter{ 0.0f };
        int edgeCount = 0;

        auto calculateFaceCentroid = [&faceCenter, &edgeCount](HalfEdgeEdge* edge)
        {
            math::vec3 pos = edge->edgePosition;
            faceCenter += pos;
            edgeCount++;
        };
        forEachEdge(calculateFaceCentroid);

        centroid = faceCenter / static_cast<float>(edgeCount);

        math::vec3& planeCentroid = centroid;
        math::vec3& planeNormal = normal;

        auto centerVerticesToCentroid = [&planeCentroid, &planeNormal](HalfEdgeEdge* edge)
        {
            float distanceToPlane = PhysicsStatics::PointDistanceToPlane(planeNormal, planeCentroid, edge->edgePosition);
            edge->edgePosition += distanceToPlane * planeNormal;
        };
        forEachEdge(calculateFaceCentroid);

        int currentEdgeId = 0;

        auto initializeEdgeToFaceFunc = [this, &currentEdgeId, edgeCount](HalfEdgeEdge* edge)
        {
            edge->face = this;

            int nextID = currentEdgeId + 1 < edgeCount ? currentEdgeId + 1 : 0;

            EdgeLabel label
            (std::make_pair(faceCount, currentEdgeId), std::make_pair(faceCount, nextID));

            edge->label = std::move(label);

            currentEdgeId++;
        };

        forEachEdge(initializeEdgeToFaceFunc);

        faceCount++;
    }

    float HalfEdgeFace::CalculateFaceArea() 
    {
        static bool notfirstTime = true;
        float sum = 0;
        math::vec3& faceCentroidPos = centroid;

        auto caculateEdgeArea = [&faceCentroidPos,&sum](HalfEdgeEdge* edge)
        {
            const math::vec3& currentEdgePos = edge->edgePosition;

            //get vector of current and next position
            math::vec3 direction = edge->nextEdge->edgePosition - currentEdgePos;
            math::vec3 directionNorm = math::normalize(direction);

            //get length of direction
            float dirVeclength = math::length(direction);
      
            //dot (centroid - current) with normalized direction
            float interpolant = math::dot(faceCentroidPos - currentEdgePos, directionNorm);

            //divide by length,this is the interpolant of the vector
            interpolant /= dirVeclength;

            //use interpolant to get closest point to centroid
            math::vec3 closestPointToCentroid = currentEdgePos + direction * interpolant;

            float triangleHeight = math::length(closestPointToCentroid - faceCentroidPos);
            //multiply length of direction with length of centroid-interpolant point divide by 2
            sum += (triangleHeight * dirVeclength * 0.5f);

            if (notfirstTime)
            {
                notfirstTime = false;
                debug::drawLine(currentEdgePos, edge->nextEdge->edgePosition, math::colors::black, 5.0f, FLT_MAX, true);
                debug::drawLine(closestPointToCentroid, faceCentroidPos, math::colors::red, 5.0f, FLT_MAX, true);
            }
        };

        forEachEdge(caculateEdgeArea);

        return sum;
    }

    void HalfEdgeFace::deleteEdges()
    {
        HalfEdgeEdge* current = startEdge->nextEdge;
        do
        {
            if (current->prevEdge->pairingEdge && current->prevEdge->pairingEdge->pairingEdge == current) current->prevEdge->pairingEdge->pairingEdge = nullptr;
            delete current->prevEdge;
            current = current->nextEdge;
        } while (current != startEdge && current != nullptr);
        startEdge = nullptr;
    }

    void HalfEdgeFace::setFaceForAllEdges()
    {
        HalfEdgeEdge* current = startEdge;
        do
        {
            current->face = this;
            current = current->nextEdge;
        } while (current != startEdge);
    }


    void HalfEdgeFace::forEachEdge(legion::core::delegate< void(HalfEdgeEdge*)> functionToExecute,
        legion::core::delegate <HalfEdgeEdge* (HalfEdgeEdge*)> getNextEdge 
    )
    {
        HalfEdgeEdge* initialEdge = startEdge;
        HalfEdgeEdge* currentEdge = startEdge;

        if (!currentEdge) return;

        //the HalfEdgeEdge* 'startEdge' creates a ring buffer.
        //This means that initialEdge will eventually go back to "startEdge", ending the loop.
        do
        {
            HalfEdgeEdge* edgeToExecuteOn = currentEdge;
            currentEdge = getNextEdge(currentEdge);
            functionToExecute(edgeToExecuteOn);

        } while (initialEdge != currentEdge && getNextEdge(currentEdge) != nullptr);

    }

    void HalfEdgeFace::forEachEdgeReverse(legion::core::delegate<void(HalfEdgeEdge*)> functionToExecute)
    {
        auto getPrevEdges = [](HalfEdgeEdge* current) {return current->prevEdge; };

        forEachEdge(functionToExecute, getPrevEdges);
    }

    void HalfEdgeFace::inverse()
    {
        normal = -normal; // Inverse the normal

        //collect edges into std::vector
        std::vector<HalfEdgeEdge*> edges;

        auto collectEdges = [&edges](HalfEdgeEdge* edge)
        {
            edges.push_back(edge);
        };

        forEachEdgeReverse(collectEdges);

        for (size_t i = 0; i < edges.size(); i++)
        {
            int nextIndex =( i + 1) % edges.size();
            int prevIndex = (i - 1) == -1 ? edges.size()-1 : (i - 1);

            HalfEdgeEdge* newNext = edges.at(prevIndex);
            HalfEdgeEdge* newPrev = edges.at(nextIndex);

            edges.at(i)->setNextAndPrevEdge(newNext, newPrev);

        }

        startEdge = edges.at(0);
    }

    bool HalfEdgeFace::testConvexity(const HalfEdgeFace& other) const
    {
        if (other == *this)
        {
            log::warn("Testing face with itself for convexity: returning true");
            return true;
        }
        math::vec3 difference = startEdge->edgePosition - other.centroid;
        float scaledAngle = math::dot(difference, normal);

        // if the scaledAngle is smaller or equal to 0, it is not convex
        if (scaledAngle <= 0)
        {
            return false;
        }

        return true;
    }

    bool HalfEdgeFace::makeNormalsConvexWithFace(HalfEdgeFace& other)
    {
        if (other == *this)
        {
            log::warn("Make normals for face convex with itself: returning false");
            return false;
        }

        face_angle_relation relation= other.getAngleRelation(*this);

        if (relation == face_angle_relation::concave)
        {
            inverse();
            return true;
        }

        if (relation == face_angle_relation::coplanar)
        {
            log::warn("COPLANAR");
            this->normal = other.normal;
            return true;

        }

        //math::vec3 difference = startEdge->edgePosition - other.centroid;
        //float scaledAngle = math::dot(difference, normal);

        //// if the scaledAngle is smaller or equal to 0, it is not convex
        //if (scaledAngle < -math::epsilon<float>())
        //{
        //    inverse();
        //    return true;
        //}
        //else if (scaledAngle < math::epsilon<float>())
        //{
        //    this->normal = other.normal;
        //    log::warn("COPLANAR");
        //    return true;
        //}
        return false;
    }

    bool HalfEdgeFace::makeNormalsConvexWithNeighbors(HalfEdgeFace& other)
    {
        HalfEdgeEdge* start = this->startEdge;
        HalfEdgeEdge* end = nullptr;

       //do
       // {
       //     log::debug(":makeNormalsConvexWithNeighbors 1 ");
       //     if (start == other.startEdge->pairingEdge)
       //     {
       //         end = start;
       //         break;
       //     }

       //     start = start->nextEdge;
       // }
       // while (start != this->startEdge);
        end = HalfEdgeFace::findMiddleEdge(*this, other);

        if (end == nullptr)
        {
            //log::debug("end == nullptr ");
            return false;
        }

        start = end;
       // end = this->startEdge;

        do
        {
            //log::debug(":makeNormalsConvexWithNeighbors 2 ");

            HalfEdgeFace* face = start->pairingEdge->face;
            face_angle_relation relation = face->getAngleRelation(*this);

            if (relation == face_angle_relation::concave)
            {
                inverse();
                //log::debug("inversing normals");
                return true;
            }
            else if (relation == face_angle_relation::convex)
            {
                return false;
            }

            //if (relation == face_angle_relation::coplanar)
            //{
            //   // log::debug("COPLANAR");
            //    this->normal = other.normal;
            //    return true;

            //}
            start = start->nextEdge;
        } while (start != end);

        //log::debug("mission failed");
        return false;
    }

    HalfEdgeFace::face_angle_relation HalfEdgeFace::getAngleRelation(const HalfEdgeFace& other)
    {
        if (other == *this)
        {
            log::warn("Calculating face angle relation between the same face, returning coplanar");
            return face_angle_relation::coplanar;
        }
  
        //log::debug("IN FUNCTION face i {} face j {}", math::to_string(normal), math::to_string(other.normal));
        float distToPlane = math::pointToPlane(centroid, other.centroid, other.normal);
        //float distToPlane2 = math::pointToPlane(other.centroid, centroid, normal);
        /*log::debug("centroid {} other.centroid {} ", centroid, other.centroid);
        log::debug("distToPlane {} ", distToPlane);
        log::debug("distToPlane2 {} ", distToPlane2);
        log::debug("angle {} ", math::angleToPlane(centroid, other.centroid, other.normal));
        log::debug("angle2 {} ", math::angleToPlane(other.centroid, centroid, normal));*/
        //distToPlane2 = distToPlane;
        //log::debug("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\");

        static float error = math::sqrt(math::epsilon<float>()) ;
       
        // if the distance to the face is negative, is it under the other face, therefore convex
        if (distToPlane <= -error)
        {
            return face_angle_relation::convex;
        }
        else if (distToPlane >= error)
        {
            return face_angle_relation::concave;
        }
        else return face_angle_relation::coplanar;
    }

    bool HalfEdgeFace::testConvexity(const HalfEdgeFace& first, const HalfEdgeFace& second)
    {
        if (first == second)
        {
            log::warn("Testing face with itself for convexity: returning true");
            return true;
        }
        return first.testConvexity(second) && second.testConvexity(first);
    }

    bool HalfEdgeFace::makeNormalsConvexWithFace(HalfEdgeFace& first, HalfEdgeFace& second)
    {
        if (first == second)
        {
            log::warn("Make normals for face convex with itself: returning false");
            return false;
        }
        bool inversedNormal = first.makeNormalsConvexWithFace(second);
        inversedNormal |= second.makeNormalsConvexWithFace(first);
        return inversedNormal;
    }

     void HalfEdgeFace::mergeCoplanarNeighbors(std::vector<HalfEdgeFace*> & removed)
     {

        HalfEdgeEdge* current = startEdge;
        if (this->normal != math::vec3(0, 0, 0))
        {
            return;
        }

        do
        {
            HalfEdgeFace* neigbor = current->pairingEdge->face;

            if (neigbor->normal != math::vec3(0, 0, 0))
            {
                current = current->nextEdge;
                continue;
            }
            face_angle_relation relation = this->getAngleRelation(*neigbor);

            if (relation == face_angle_relation::coplanar)
            {
                //log::debug("relation == face_angle_relation::coplanar in :mergeCoplanarNeighbors ");
                HalfEdgeEdge* middle = HalfEdgeFace::findMiddleEdge(*this, *neigbor);
                if (middle)
                {
                    this->mergeFaces(*middle);

                    removed.push_back(neigbor);
                    neigbor->normal = math::vec3();
                }
           
            }

           // this->mergeFaces
            current = current->nextEdge;
        }
        while (current != startEdge);

        return;
    }



    HalfEdgeEdge* HalfEdgeFace::findMiddleEdge(const HalfEdgeFace& first, const HalfEdgeFace& second)
    {
        HalfEdgeEdge* firstCurrent = first.startEdge;
        do
        {
            HalfEdgeEdge* secondCurrent = second.startEdge;
            do
            {
                if (firstCurrent->pairingEdge == secondCurrent)
                {
                    return firstCurrent;
                }
                if (firstCurrent == secondCurrent->pairingEdge)
                {
                    return firstCurrent;
                }


                secondCurrent = secondCurrent->nextEdge;
            } while (secondCurrent != second.startEdge);
            firstCurrent = firstCurrent->nextEdge;
        } while (firstCurrent != first.startEdge);

        return nullptr;
    }

    HalfEdgeFace* HalfEdgeFace::mergeFaces(HalfEdgeEdge& middleEdge)
    {
        if (middleEdge.face->startEdge == &middleEdge)
        {
            middleEdge.face->startEdge = middleEdge.face->startEdge->prevEdge;
        }
        // Set correct faces
        HalfEdgeEdge* end = middleEdge.pairingEdge;
        HalfEdgeEdge* current = middleEdge.pairingEdge;
        do
        {
            current->face = middleEdge.face;
            current = current->nextEdge;
        } while (current != end);

        // Link edges
        middleEdge.prevEdge->nextEdge = middleEdge.pairingEdge->nextEdge;
        middleEdge.nextEdge->prevEdge = middleEdge.pairingEdge->prevEdge;
        middleEdge.pairingEdge->prevEdge->nextEdge = middleEdge.nextEdge;
        middleEdge.pairingEdge->nextEdge->prevEdge = middleEdge.prevEdge;

        HalfEdgeFace* face = middleEdge.face;

        // Recalculate centroid
        {
            math::vec3 faceCenter{ 0.0f };
            int edgeCount = 0;

            auto calculateFaceCentroid = [&faceCenter, &edgeCount](HalfEdgeEdge* edge)
            {
                math::vec3 pos = edge->edgePosition;
                faceCenter += pos;
                edgeCount++;
            };
            face->forEachEdge(calculateFaceCentroid);

            face->centroid = faceCenter / static_cast<float>(edgeCount);
        }

        delete middleEdge.pairingEdge;
        delete& middleEdge;

        return face;
    }

    void HalfEdgeFace::DEBUG_DrawFace(const math::mat4& transform,const math::color& debugColor,float time)
    {

        auto drawFunc = [&transform,debugColor,time](HalfEdgeEdge* edge)
        {
            edge->DEBUG_drawEdge(transform, debugColor, time);

        };

        forEachEdge(drawFunc);
    }

    void HalfEdgeFace::DEBUG_DirectionDrawFace(const math::mat4& transform, const math::color& debugColor, float time)
    {
        auto drawFunc = [&transform, debugColor, time](HalfEdgeEdge* edge)
        {
            edge->DEBUG_directionDrawEdge(transform, debugColor, time, 5.0f);
        };

        forEachEdge(drawFunc);
    }

    HalfEdgeFace::~HalfEdgeFace()
    {
        auto deleteFunc = [](HalfEdgeEdge* edge)
        {
            if (edge->pairingEdge && edge->pairingEdge && edge->pairingEdge->pairingEdge == edge)
            {
                edge->pairingEdge->pairingEdge = nullptr;
            }
            delete edge;
        };

        forEachEdge(deleteFunc);

    }
}
