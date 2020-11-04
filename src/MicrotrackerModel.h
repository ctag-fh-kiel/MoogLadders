// Based on an implementation by Magnus Jonsson
// https://github.com/magnusjonsson/microtracker (unlicense)

#pragma once

#ifndef MICROTRACKER_MODEL_H
#define MICROTRACKER_MODEL_H

#include "LadderFilterBase.h"
#include "util.h"
#include "helpers/ctagFastMath.hpp"

namespace Moog {
    class MicrotrackerMoog : public LadderFilterBase {

    public:

        MicrotrackerMoog(float sampleRate) : LadderFilterBase(sampleRate) {
            p0 = p1 = p2 = p3 = p32 = p33 = p34 = 0.0;
            SetCutoff(1000.0f);
            SetResonance(0.10f);
        }

        virtual ~MicrotrackerMoog() {}

        virtual void Process(float *samples, uint32_t n) override {
            float k = resonance * 4;
            for (int s = 0; s < n; ++s) {
                // Coefficients optimized using differential evolution
                // to make feedback gain 4.0 correspond closely to the
                // border of instability, for all values of omega.
                float out = p3 * 0.360891 + p32 * 0.417290 + p33 * 0.177896 + p34 * 0.0439725;

                p34 = p33;
                p33 = p32;
                p32 = p3;

                float tanhp0 = CTAG::SP::HELPERS::fasttanh(p0);
                p0 += (fast_tanh(samples[s] - k * out) - tanhp0) * cutoff;
                float tanhp1 = CTAG::SP::HELPERS::fasttanh(p1);
                p1 += (tanhp0 - tanhp1) * cutoff;
                float tanhp2 = CTAG::SP::HELPERS::fasttanh(p2);
                p2 += (tanhp1 - tanhp2) * cutoff;
                p3 += (tanhp2 - fast_tanh(p3)) * cutoff;

                samples[s] = out;
            }
        }

        virtual void SetResonance(float r) override {
            resonance = r;
        }

        virtual void SetCutoff(float c) override {
            cutoff = c * 2 * MOOG_PI / sampleRate;
            cutoff = moog_min(cutoff, 1);
        }

    private:

        float p0;
        float p1;
        float p2;
        float p3;
        float p32;
        float p33;
        float p34;
    };
}
#endif
