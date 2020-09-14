// See: http://www.willpirkle.com/forum/licensing-and-book-code/licensing-and-using-book-code/
// The license is "You may also use the code from the FX and Synth books without licensing or fees. 
// The code is for you to develop your own plugins for your own use or for commercial use."

#pragma once

#ifndef OBERHEIM_VARIATION_LADDER_H
#define OBERHEIM_VARIATION_LADDER_H

#include "LadderFilterBase.h"
#include "util.h"
#include "helpers/ctagFastMath.hpp"

class VAOnePole
{
public:

	VAOnePole(float sr) : sampleRate(sr)
	{
		Reset();
	}
	
	void Reset()
	{
		alpha = 1.0;
		beta = 0.0;
		gamma = 1.0;
		delta = 0.0;
		epsilon = 0.0;
		a0 = 1.0;
		feedback = 0.0;
		z1 = 0.0;
	}
	
	float Tick(float s)
	{
		s = s * gamma + feedback + epsilon * GetFeedbackOutput();
		float vn = (a0 * s - z1) * alpha;
		float out = vn + z1;
		z1 = vn + out;
		return out;
	}
	
	void SetFeedback(float fb) { feedback = fb; }
	float GetFeedbackOutput(){ return beta * (z1 + feedback * delta); }
	void SetAlpha(float a) { alpha = a; };
	void SetBeta(float b) { beta = b; };
	
private:

	float sampleRate;
	float alpha;
	float beta;
	float gamma;
	float delta;
	float epsilon;
	float a0;
	float feedback;
	float z1;
};

class OberheimVariationMoog : public LadderFilterBase
{
	
public:
	
	OberheimVariationMoog(float sampleRate) : LadderFilterBase(sampleRate)
	{
		LPF1 = new VAOnePole(sampleRate);
		LPF2 = new VAOnePole(sampleRate);
		LPF3 = new VAOnePole(sampleRate);
		LPF4 = new VAOnePole(sampleRate);
		
		saturation = 1.0;
		Q = 3.0;
		
		SetCutoff(1000.f);
		SetResonance(0.1f);
	}
	
	virtual ~OberheimVariationMoog()
	{
		delete LPF1;
		delete LPF2;
		delete LPF3;
		delete LPF4;
	}
	
	virtual void Process(float * samples, uint32_t n) noexcept override
	{
		for (int s = 0; s < n; ++s)
		{
			float input = samples[s];
			
			float sigma =
				LPF1->GetFeedbackOutput() +
				LPF2->GetFeedbackOutput() +
				LPF3->GetFeedbackOutput() +
				LPF4->GetFeedbackOutput();
			
			input *= 1.0 + K;
			
			// calculate input to first filter
			float u = (input - K * sigma) * alpha0;
			
			u = tanh(saturation * u);
			
			float stage1 = LPF1->Tick(u);
			float stage2 = LPF2->Tick(stage1);
			float stage3 = LPF3->Tick(stage2);
			float stage4 = LPF4->Tick(stage3);
			
			// Oberheim variations
			samples[s] =
				oberheimCoefs[0] * u +
				oberheimCoefs[1] * stage1 +
				oberheimCoefs[2] * stage2 +
				oberheimCoefs[3] * stage3 +
				oberheimCoefs[4] * stage4;
		}
	}
	
	virtual void SetResonance(float r) override
        {
             // this maps resonance = 1->10 to K = 0 -> 4
             K = (4.0) * (r - 1.0)/(10.0 - 1.0);
        }

	virtual void SetCutoff(float c) override
	{
		cutoff = c;
		
		// prewarp for BZT
		float wd = 2.0 * MOOG_PI * cutoff;
		float T = 1.0 / sampleRate;
		float wa = (2.0 / T) * CTAG::SP::HELPERS::fasttan(wd * T / 2.0);
		float g = wa * T / 2.0;
		
		// Feedforward coeff
		float G = g / (1.0 + g);
		
		LPF1->SetAlpha(G);
		LPF2->SetAlpha(G);
		LPF3->SetAlpha(G);
		LPF4->SetAlpha(G);

		LPF1->SetBeta(G*G*G / (1.0 + g));
		LPF2->SetBeta(G*G / (1.0 + g));
		LPF3->SetBeta(G / (1.0 + g));
		LPF4->SetBeta(1.0 / (1.0 + g));
		
		gamma = G*G*G*G;
		alpha0 = 1.0 / (1.0 + K * gamma);
		
		// Oberheim variations / LPF4
		oberheimCoefs[0] = 0.0;
		oberheimCoefs[1] = 0.0;
		oberheimCoefs[2] = 0.0;
		oberheimCoefs[3] = 0.0;
		oberheimCoefs[4] = 1.0;
	}
	
private:
	
	VAOnePole * LPF1;
	VAOnePole * LPF2;
	VAOnePole * LPF3;
	VAOnePole * LPF4;
	
	float K;
	float gamma;
	float alpha0;
	float Q;
	float saturation;
	
	float oberheimCoefs[5];
};

#endif
