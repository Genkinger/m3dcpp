#include "m3dcpp.h"

#include <iostream>
namespace M3D {
    bool InRange(float value, float min, float max) {
        return (value >= min && value <= max);
    }

    float Map(float value, float in_min, float in_max, float out_min, float out_max) {
        float val = Constrain(value,in_min,in_max);

        return (val - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    }

    float Constrain(float value, float min, float max) {
        if (value <= max && value >= min) {
            return value;
        }else if (value < min) {
            return min;
        }else if (value > max) {
            return max;
        }
    }

    float InterpolateLinear(float start, float end, float factor) {
        return (1 - factor) * start + factor * end;
    }
}