#include <physics/halfedgeedge.hpp>
#include <physics/halfedgeface.hpp>

namespace legion::physics
{
    HalfEdgeEdge::HalfEdgeEdge(math::vec3 newEdgePositionPtr) : edgePosition{ newEdgePositionPtr }
    {
        static int idCount = 0;
        identifier = idCount;

        if (identifier == 433)
        {
            //DebugBreak();
        }
        idCount++;
    }
    void HalfEdgeEdge::setNextAndPrevEdge(HalfEdgeEdge* newPrevEdge, HalfEdgeEdge* newNextEdge)
    {
        nextEdge = newNextEdge;
        prevEdge = newPrevEdge;
    }

    void HalfEdgeEdge::setNext(HalfEdgeEdge* newNextEdge)
    {
        nextEdge = newNextEdge;
        newNextEdge->prevEdge = this;

    }

    void HalfEdgeEdge::setPrev(HalfEdgeEdge* newPrevEdge)
    {
        prevEdge = newPrevEdge;
        newPrevEdge->nextEdge = this;
    }

    void HalfEdgeEdge::setPairingEdge(HalfEdgeEdge* edge)
    {
        pairingEdge = edge;
        edge->pairingEdge = this;
    }

    math::vec3 HalfEdgeEdge::getLocalNormal() const
    {
        return face->normal;
    }

    math::vec3 HalfEdgeEdge::getLocalEdgeDirection() const
    {
        return nextEdge->edgePosition - edgePosition;
    }

    bool HalfEdgeEdge::isVertexVisible(const math::vec3& vert,float epsilon)
    {
        float distanceToPlane =
            math::pointToPlane(vert, face->centroid, face->normal);

        return distanceToPlane > epsilon;
    }

    bool HalfEdgeEdge::isEdgeHorizonFromVertex(const math::vec3& vert, float epsilon)
    {
        return isVertexVisible(vert,epsilon) && !pairingEdge->isVertexVisible(vert,epsilon);
    }

    void HalfEdgeEdge::DEBUG_drawEdge(const math::mat4& transform, const math::color& debugColor, float time, float width)
    {
        math::vec3 worldStart = transform * math::vec4(edgePosition, 1);
        math::vec3 worldEnd = transform * math::vec4(nextEdge->edgePosition, 1);

        debug::user_projectDrawLine(worldStart, worldEnd, debugColor, width, time, true);
    }

    void HalfEdgeEdge::DEBUG_drawInsetEdge(const math::vec3 spacing, const math::color& debugColor, float time, float width)
    {
        math::vec3 worldCentroid = face->centroid + spacing;

        math::vec3 worldStart = edgePosition + spacing;
        math::vec3 startDifference = (worldCentroid - worldStart) * 0.1f;

        math::vec3 worldEnd = nextEdge->edgePosition + spacing;
        math::vec3 endDifference = (worldCentroid - worldEnd) * 0.1f;

        debug::user_projectDrawLine(worldStart + startDifference, worldEnd + endDifference, debugColor, width, time, true);
    }

    void HalfEdgeEdge::DEBUG_directionDrawEdge(const math::mat4& transform, const math::color& debugColor, float time, float width)
    {
        math::vec3 worldStart = transform * math::vec4(edgePosition, 1);
        math::vec3 worldEnd = transform * math::vec4(nextEdge->edgePosition, 1);

        math::vec3 worldCentroid = transform * math::vec4(face->centroid, 1);

        math::vec3 startDifference = (worldCentroid - worldStart) * 0.1f;
        math::vec3 endDifference = (worldCentroid - worldEnd) * 0.1f;

        debug::user_projectDrawLine(worldStart + startDifference, worldEnd + endDifference, debugColor, width, time, true);

        math::vec3 pointStart = worldStart + startDifference;
        math::vec3 diff = worldEnd + endDifference - (worldStart + startDifference);

        debug::user_projectDrawLine(pointStart + diff * 0.75f, pointStart + diff * 0.75f + math::vec3(-0.1f,0,0), math::colors::red, width, time, true);
    }

    void HalfEdgeEdge::suicidalMergeWithPairing(math::mat4 transform,  bool shouldDebug)
    {
        //----//
        HalfEdgeFace* mergeFace = pairingEdge->face;
        HalfEdgeEdge* prevFromCurrent = pairingEdge->prevEdge->pairingEdge->face == this->face ? pairingEdge->nextEdge : pairingEdge->prevEdge;

        //next-pairing->next
        HalfEdgeEdge* prevFromCurrentConnection = prevFromCurrent->nextEdge->pairingEdge->nextEdge;
        prevFromCurrent->setNext(prevFromCurrentConnection);

        HalfEdgeEdge* nextFromCurrent = pairingEdge->nextEdge->pairingEdge->face == this->face ? pairingEdge->prevEdge : pairingEdge->nextEdge;
        HalfEdgeEdge* nextFromCurrentConnection = nextFromCurrent->prevEdge->pairingEdge->prevEdge;

        nextFromCurrent->setPrev(nextFromCurrentConnection);

        if (shouldDebug)
        {
            prevFromCurrent->DEBUG_drawEdge(transform, math::colors::red,FLT_MAX);
            prevFromCurrentConnection->DEBUG_drawEdge(transform, math::colors::blue, FLT_MAX);

            nextFromCurrent->DEBUG_drawEdge(transform, math::colors::green, FLT_MAX);
            nextFromCurrentConnection->DEBUG_drawEdge(transform, math::colors::magenta, FLT_MAX);

        }

        face->startEdge = prevFromCurrent;
        face->initializeFace();

        if (prevEdge->pairingEdge->face == mergeFace)
        {
            delete prevEdge->pairingEdge;
            delete prevEdge;
        }

        if (nextEdge->pairingEdge->face == mergeFace)
        {
            delete nextEdge->pairingEdge;
            delete nextEdge;
        }

        mergeFace->startEdge = nullptr;
        delete mergeFace;

        delete pairingEdge;
        delete this;

    }
}

