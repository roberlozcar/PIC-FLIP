#ifndef PARTICLES_2_H
#define PARTICLES_2_H

#include "Math/Vector2.h"
#include "Math/Vector3.h"
#include <vector>

namespace asa
{
class Particles2
{
public:
    Particles2(const uint n = 0)
    {
        size = n;
        position.resize(size);
        velocity.resize(size);
        ink.resize(size);
    }

    uint getSize() const { return size; }

    void resize(const uint n)
    {
        if (size > n) {
            size = n;
            position.resize(size);
            velocity.resize(size);
            ink.resize(size);
        } else if (size < n) {
            position.resize(n);
            velocity.resize(n);
            ink.resize(n);

            const uint diff = n - size;
            std::memset(&position[size], 0, sizeof(Vector2) * diff);
            std::memset(&velocity[size], 0, sizeof(Vector2) * diff);
            std::memset(&ink[size], 0, sizeof(Vector3) * diff);
            size = n;
        }
    }

    const Vector2 &getPosition(const uint i) const { return position[i]; }

    const Vector2 &getVelocity(const uint i) const { return velocity[i]; }

    const Vector3 &getInk(const uint i) const { return ink[i]; }

    void setPosition(const uint i, const Vector2 &pos) { position[i] = pos; }

    void setVelocity(const uint i, const Vector2 &vel) { velocity[i] = vel; }

    void setInk(const uint i, const Vector3 s) { ink[i] = s; }

    uint addParticle(const Vector2 &pos, const Vector2 &vel = Vector2(), const Vector3 &s = Vector3())
    {
        const uint id = size++;

        position.push_back(pos);
        velocity.push_back(vel);
        ink.push_back(s);

        return id;
    }

private:
    uint size;
    std::vector<Vector2> position;
    std::vector<Vector2> velocity;
    std::vector<Vector3> ink;
};
}  // namespace asa

#endif
