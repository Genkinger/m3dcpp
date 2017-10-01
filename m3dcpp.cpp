#include "m3dcpp.h"

#include <iostream>
namespace M3D {
    bool InRange(float value, float min, float max) {
        if (value >= min && value <= max) {
            return true;
        } else {
            return false;
        }
    }
    float Map(float value, float in_min, float in_max, float out_min, float out_max) {
        return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    }
    float Constrain(float value, float min, float max) {
        if (InRange(value, min, max)) {
            return value;
        }else if (value <= min) {
            return min;
        }else if (value >= max) {
            return max;
        }
        return value; // get rid of compiler warning
    }
    float InterpolateLinear(float start, float end, float factor) {
        return (1 - factor) * start + factor * end;
    }
}